#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Basis.h"
#include "Grid.h"
#include "Field.h"
#include "CFA.h"
#include "LinAlg.h"
#include "CDG.h"

#define ADV_FORWARD (+1)
#define ADV_BACKWARD (-1)

//#define CDG_TEST 1

CDG::CDG( Field* _phi, Field* _velx, Field* _vely ) : CFA( _phi, _velx, _vely ) {
	betaInv_ij = NULL;
}

CDG::~CDG() {
	int i;

	for( i = 0; i < phi->grid->nPolys; i++ ) {
		delete[] betaInv_ij[i];
	}
	delete[] betaInv_ij;
}

void CDG::InitBetaIJInv( Func* func ) {
	int 		i, j, k, l, pi;
	Grid* 		grid 		= phi->grid;
	Polygon*	poly;
	Triangle*	tri;
	Basis*		basis;
	int			nBasis		= phi->basis[0]->nFuncs;
	double		fj[nBasis];
	double		weight, *coord;
	double		beta_ij[nBasis*nBasis];
	double		volErr;

	if( betaInv_ij == NULL ) {
		betaInv_ij = new double*[grid->nPolys];
		for( pi = 0; pi < grid->nPolys; pi++ ) {
			betaInv_ij[pi] = NULL;
		}
	}

	/* set up the matrix inverse and intial basis coefficients */
	for( pi = 0; pi < grid->nPolys; pi++ ) {
		if( betaInv_ij[pi] == NULL ) {
			betaInv_ij[pi] = new double[nBasis*nBasis];
		}

		for( j = 0; j < nBasis*nBasis; j++ ) {
			beta_ij[j] = 0.0;
			betaInv_ij[pi][j] = 0.0;
		}

		poly = grid->polys[pi];
		basis = phi->basis[pi];
		for( k = 0; k < poly->n; k++ ) {
			tri = poly->tris[k];
			for( l = 0; l < tri->nQuadPts; l++ ) {
				for( j = 0; j < nBasis; j++ ) {
					for( i = 0; i < nBasis; i++ ) {
						weight = tri->wi[l]*tri->Area()/poly->Area();
						coord = tri->qi[l];
						beta_ij[j*nBasis+i] += weight*basis->EvalIJ( coord, i )*basis->EvalIJ( coord, j );
					}
				}
			}
		}
		MatInv( beta_ij, betaInv_ij[pi], nBasis );
		for( j = 0; j < nBasis; j++ ) {
			fj[j] = 0.0;
			for( k = 0; k < poly->n; k++ ) {
				tri = poly->tris[k];
				for( l = 0; l < tri->nQuadPts; l++ ) {
					weight = tri->wi[l]*tri->Area()/poly->Area();
					coord = tri->qi[l];
					/* basis initially set as the spatial values at the poly coordinates */
					fj[j] += weight*basis->EvalIJ( coord, j )*func( coord );
				}
			}
		}
		/* set the initial basis coefficients */
		AXEB( betaInv_ij[pi], fj, basis->ci, nBasis );

		if( !basis->TestMean( &volErr ) ) {
			cout << "ERROR: basis function mean not equal to first component..." << volErr << endl;
		}
	}
}

void CDG::Advect( double dt ) {
	Grid* 	grid	= phi->grid;
	Grid*	preGrid = new Grid( grid->nx, grid->ny, grid->minx, grid->miny, grid->maxx, grid->maxy, grid->quadOrder, grid->basisOrder, true );
	Field*	phiTemp = new Field( grid );

	CalcChars( preGrid, dt );
	preGrid->UpdateEdges();
	preGrid->UpdatePolys();
	preGrid->UpdateTriangles();
	CalcFluxes( preGrid, phiTemp, dt );
	phi->Copy( phiTemp );

	delete phiTemp;
	delete preGrid;
}

void CDG::BasisProjection( int kp, int k, double* Pij ) {
	int tri_i, quad_i, basis_m, basis_j;
	Grid* 		grid 		= phi->grid;
	Polygon* 	poly 		= grid->polys[kp];
	Basis*		basis_k		= phi->basis[k];
	Basis*		basis_kp	= phi->basis[kp];
	Triangle* 	tri;
	int			nBasis 		= grid->basisOrder*grid->basisOrder;
	int			nBasis2		= nBasis*nBasis;
	double		beta_mj[nBasis2];
	double		weight, *coord;

	for( basis_j = 0; basis_j < nBasis2; basis_j++ ) {
		beta_mj[basis_j] = 0.0;
	}

	for( tri_i = 0; tri_i < poly->n; tri_i++ ) {
		tri = poly->tris[tri_i];
		for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
			weight = tri->wi[quad_i]*tri->Area()/poly->Area();
			coord = tri->qi[quad_i];
			for( basis_m = 0; basis_m < nBasis; basis_m++ ) {
				for( basis_j = 0; basis_j < nBasis; basis_j++ ) {
					beta_mj[basis_m*nBasis+basis_j] += weight*basis_kp->EvalIJ( coord, basis_m )*basis_k->EvalIJ( coord, basis_j );
				}
			}
		}
	}

	Mult( betaInv_ij[kp], beta_mj, Pij, nBasis );
}

void CDG::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 		poly_i, edge_i, basis_i, tri_i, quad_i;
	Grid*		grid	= phi->grid;
	Polygon 	*prePoly, *intPoly, *incPoly;
	Triangle*	tri;
	int			pinds[6], into, from;
	double 		weight, tracer, basis_into, basis_from;
	double**	flux, qf[2];
	int			nBasis	= phi->basis[0]->nFuncs;

	flux = new double*[grid->nPolys];
	for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
		flux[poly_i] = new double[nBasis];
		for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
			flux[poly_i][basis_i] = 0.0;
		}
	}

	for( edge_i = 0; edge_i < grid->nEdges; edge_i++ ) {
        prePoly = CreatePreImage( edge_i, grid, preGrid, &into, &from, pinds );
		if( prePoly == NULL ) {
			continue;
		}

        for( poly_i = 0; poly_i < 6; poly_i++ ) {
			from = pinds[poly_i];
			incPoly = grid->polys[from];
			intPoly = Intersection( prePoly, incPoly );

			if( intPoly ) {
				for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
					tri = intPoly->tris[tri_i];
					for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
						TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
						weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
						tracer = phi->EvalAtCoord( tri->qi[quad_i] );
						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							basis_into = phi->basis[into]->EvalIJ( qf, basis_i );
							basis_from = phi->basis[from]->EvalIJ( qf, basis_i );
							flux[into][basis_i] += weight*tracer*basis_into;
							flux[from][basis_i] -= weight*tracer*basis_from;
						}
					}
				}
				delete intPoly;
			}
        }
        delete prePoly;
	}

	/* add rhs contributions from previous poly */
	for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
		for( tri_i = 0; tri_i < grid->polys[poly_i]->n; tri_i++ ) {
			tri = grid->polys[poly_i]->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
				weight = tri->wi[quad_i]*tri->Area()/grid->polys[poly_i]->Area();
				tracer = phi->EvalAtCoord( tri->qi[quad_i] );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					basis_into = phi->basis[poly_i]->EvalIJ( qf, basis_i );
					flux[poly_i][basis_i] += weight*tracer*basis_into;
				}
			}
		}

		/* update the poly coefficients */
		AXEB( betaInv_ij[poly_i], flux[poly_i], phiTemp->basis[poly_i]->ci, nBasis );
	}

	if( phiTemp->basis[0]->order == 2 ) {
		;//Limiter_FirstOrder( phiTemp );
	}
	if( phiTemp->basis[0]->order > 2 ) {
		Limiter_SecondOrder( phiTemp );
	}

	for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
		delete[] flux[poly_i];
	}
	delete[] flux;
}

void CDG::Limiter_FirstOrder( Field* phiTemp ) {
	Grid*	grid	= phi->grid;
	Basis*	basis;
	int		nx, ny;
	int		poly_i, basis_i;
	double	alpha, dxTmp, dyTmp;

	for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
		basis = phiTemp->basis[poly_i];
		nx = poly_i%grid->nx;
		ny = poly_i/grid->nx;

		if( nx > 0 && nx < grid->nx - 1 && ny > 0 && ny < grid->ny - 1 ) {
			alpha = PolyLimiter( phiTemp, poly_i );

			if( fabs( alpha ) > 1.0 - 1.0e-4 ) {
				continue;
			}

			dxTmp = basis->ci[1];
			dyTmp = basis->ci[basis->order];

			for( basis_i = 1; basis_i < basis->nFuncs; basis_i++ ) {
				basis->ci[basis_i] = 0.0;
			}
			basis->ci[1]            = alpha*dxTmp;
			basis->ci[basis->order] = alpha*dyTmp;
		}
	}
}

void CDG::Limiter_SecondOrder( Field* phiTemp ) {
	Grid*	grid	= phi->grid;
	Basis*	basis;
	int		nx, ny;
	int		poly_i, basis_i;
	double	alpha_1, alpha_x, alpha_y, alpha_2, dxTmp, dyTmp, dxxTmp, dyyTmp, dxyTmp;

	for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
		basis = phiTemp->basis[poly_i];
		nx = poly_i%grid->nx;
		ny = poly_i/grid->nx;

		if( nx > 0 && nx < grid->nx - 1 && ny > 0 && ny < grid->ny - 1 ) {
			alpha_1 = PolyLimiter( phiTemp, poly_i );
			alpha_x = PolySecondLimiter( phiTemp, poly_i, 0 );
			alpha_y = PolySecondLimiter( phiTemp, poly_i, 1 );

			alpha_2 = ( alpha_x < alpha_y ) ? alpha_x : alpha_y;
			alpha_1 = ( alpha_1 > alpha_2 ) ? alpha_1 : alpha_2;

			if( fabs( alpha_1 ) > 1.0 - 1.0e-4 && fabs( alpha_2 ) > 1.0 - 1.0e-4 ) {
				continue;
			}

			dxTmp = basis->ci[1];
			dyTmp = basis->ci[basis->order];
			dxxTmp = basis->ci[2];
			dyyTmp = basis->ci[2*basis->order];
			dxyTmp = basis->ci[basis->order+1];

			for( basis_i = 1; basis_i < basis->nFuncs; basis_i++ ) {
				basis->ci[basis_i] = 0.0;
			}
			basis->ci[1]              = alpha_1*dxTmp;
			basis->ci[basis->order]   = alpha_1*dyTmp;
			basis->ci[2]              = alpha_2*dxxTmp;
			basis->ci[2*basis->order] = alpha_2*dyyTmp;
			basis->ci[basis->order+1] = alpha_2*dxyTmp;
		}
	}
}

double CDG::PolyLimiter( Field* phiTemp, int poly_i ) {
	Grid*		grid			= phi->grid;
	Polygon* 	poly 			= grid->polys[poly_i];
	Basis*		basis 			= phiTemp->basis[poly_i];
	int			vinds[4];
	int			nx				= poly_i%grid->nx;
	int			ny				= poly_i/grid->nx;
	int			vi, ci, cj;
	double		maxPhiAtVert	= -1.0e+99;
	double		minPhiAtVert	= +1.0e+99;
	double		maxPhiInPoly	= -1.0e+99;
	double		minPhiInPoly	= +1.0e+99;
	double		phiAtVert, phiInPoly, alpha, gamma;

	grid->GetPolyVertInds( poly_i, vinds );

	for( vi = 0; vi < poly->n; vi++ ) {
		phiAtVert = basis->EvalFull( poly->verts[vi] ) - basis->ci[0];
		minPhiAtVert = ( phiAtVert < minPhiAtVert ) ? phiAtVert : minPhiAtVert;
		maxPhiAtVert = ( phiAtVert > maxPhiAtVert ) ? phiAtVert : maxPhiAtVert;
	}
	for( cj = ny - 1; cj < ny + 2; cj++ ) {
		for( ci = nx - 1; ci < nx + 2; ci++ ) {
			if( ci == nx && cj == ny ) {
				continue;
			}
			phiInPoly = phiTemp->basis[cj*grid->nx+ci]->ci[0] - basis->ci[0];
			minPhiInPoly = ( phiInPoly < minPhiInPoly ) ? phiInPoly : minPhiInPoly;
			maxPhiInPoly = ( phiInPoly > maxPhiInPoly ) ? phiInPoly : maxPhiInPoly;
		}
	}

	alpha = 1.0;
	if( fabs( minPhiAtVert ) > 1.0e-6 ) {
		alpha = minPhiInPoly/minPhiAtVert > 0.0 ? minPhiInPoly/minPhiAtVert : 0.0;
	}

	gamma = 1.0;
	if( fabs( maxPhiAtVert ) > 1.0e-6 ) {
		gamma = maxPhiInPoly/maxPhiAtVert > 0.0 ? maxPhiInPoly/maxPhiAtVert : 0.0;
	}

	alpha = ( alpha < gamma ) ? alpha : gamma;
	alpha = ( alpha < 1.0 )   ? alpha : 1.0;

	return alpha;
}

double CDG::PolySecondLimiter( Field* phiTemp, int poly_i, int dim ) {
	Grid*		grid			= phi->grid;
	Polygon* 	poly 			= grid->polys[poly_i];
	Basis*		basis 			= phiTemp->basis[poly_i];
	int			vinds[4];
	int			nx				= poly_i%grid->nx;
	int			ny				= poly_i/grid->nx;
	int			vi, ci, cj, co;
	double		maxDPhiAtVert	= -1.0e+99;
	double		minDPhiAtVert	= +1.0e+99;
	double		maxDPhiInPoly	= -1.0e+99;
	double		minDPhiInPoly	= +1.0e+99;
	double		dPhiAtVert, dPhiInPoly;
	double 		alpha, gamma;

	grid->GetPolyVertInds( poly_i, vinds );

	for( vi = 0; vi < poly->n; vi++ ) {
		dPhiAtVert = basis->EvalDerivFull( poly->verts[vi], dim );
		minDPhiAtVert = ( dPhiAtVert < minDPhiAtVert ) ? dPhiAtVert : minDPhiAtVert;
		maxDPhiAtVert = ( dPhiAtVert > maxDPhiAtVert ) ? dPhiAtVert : maxDPhiAtVert;
	}
	for( cj = ny - 1; cj < ny + 2; cj++ ) {
		for( ci = nx - 1; ci < nx + 2; ci++ ) {
			if( ci == nx && cj == ny ) {
				continue;
			}
			co = cj*grid->nx + ci;
			dPhiInPoly = phiTemp->basis[co]->EvalDerivFull( grid->polys[co]->origin, dim );
			minDPhiInPoly = ( dPhiInPoly < minDPhiInPoly ) ? dPhiInPoly : minDPhiInPoly;
			maxDPhiInPoly = ( dPhiInPoly > maxDPhiInPoly ) ? dPhiInPoly : maxDPhiInPoly;
		}
	}

	alpha = 1.0;
	if( fabs( minDPhiAtVert ) > 1.0e-6 ) {
		alpha = minDPhiInPoly/minDPhiAtVert > 0.0 ? minDPhiInPoly/minDPhiAtVert : 0.0;
	}

	gamma = 1.0;
	if( fabs( maxDPhiAtVert ) > 1.0e-6 ) {
		gamma = maxDPhiInPoly/maxDPhiAtVert > 0.0 ? maxDPhiInPoly/maxDPhiAtVert : 0.0;
	}

	alpha = ( alpha < gamma ) ? alpha : gamma;
	alpha = ( alpha < 1.0 )   ? alpha : 1.0;

	return alpha;
}
