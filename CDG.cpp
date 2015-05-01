#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
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

	for( i = 0; i < phi->grid->nCells; i++ ) {
		delete[] betaInv_ij[i];
	}
	delete[] betaInv_ij;
}

void CDG::InitBetaIJInv( Func* func ) {
	int 		i, j, k, l, pi;
	Grid* 		grid 		= phi->grid;
	Cell*		cell;
	Triangle*	tri;
	Basis*		basis;
	int			nBasis		= phi->basis[0]->nFuncs;
	double		fj[nBasis];
	double		weight, *coord;
	double		beta_ij[nBasis*nBasis];
	double		volErr;

	if( betaInv_ij == NULL ) {
		betaInv_ij = new double*[grid->nCells];
		for( pi = 0; pi < grid->nCells; pi++ ) {
			betaInv_ij[pi] = NULL;
		}
	}

	/* set up the matrix inverse and intial basis coefficients */
	for( pi = 0; pi < grid->nCells; pi++ ) {
		if( betaInv_ij[pi] == NULL ) {
			betaInv_ij[pi] = new double[nBasis*nBasis];
		}

		for( j = 0; j < nBasis*nBasis; j++ ) {
			beta_ij[j] = 0.0;
			betaInv_ij[pi][j] = 0.0;
		}

		cell = grid->cells[pi];
		basis = phi->basis[pi];
		for( k = 0; k < cell->n; k++ ) {
			tri = cell->tris[k];
			for( l = 0; l < tri->nQuadPts; l++ ) {
				for( j = 0; j < nBasis; j++ ) {
					for( i = 0; i < nBasis; i++ ) {
						weight = tri->wi[l]*tri->Area()/cell->Area();
						coord = tri->qi[l];
						beta_ij[j*nBasis+i] += weight*basis->EvalIJ( coord, i )*basis->EvalIJ( coord, j );
					}
				}
			}
		}
		MatInv( beta_ij, betaInv_ij[pi], nBasis );
		for( j = 0; j < nBasis; j++ ) {
			fj[j] = 0.0;
			for( k = 0; k < cell->n; k++ ) {
				tri = cell->tris[k];
				for( l = 0; l < tri->nQuadPts; l++ ) {
					weight = tri->wi[l]*tri->Area()/cell->Area();
					coord = tri->qi[l];
					/* basis initially set as the spatial values at the cell coordinates */
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
	preGrid->UpdateCells();
	preGrid->UpdateTriangles();
	CalcFluxes( preGrid, phiTemp, dt );
	phi->Copy( phiTemp );

	delete phiTemp;
	delete preGrid;
}

void CDG::BasisProjection( int kp, int k, double* Pij ) {
	int tri_i, quad_i, basis_m, basis_j;
	Grid* 		grid 		= phi->grid;
	Cell* 		cell 		= grid->cells[kp];
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

	for( tri_i = 0; tri_i < cell->n; tri_i++ ) {
		tri = cell->tris[tri_i];
		for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
			weight = tri->wi[quad_i]*tri->Area()/cell->Area();
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
	int 		cell_i, edge_i, basis_i, tri_i, quad_i;
	Grid*		grid	= phi->grid;
	Polygon 	*prePoly, *intPoly;
	Cell		*incPoly;
	Triangle*	tri;
	int			pinds[6], into, from;
	double 		weight, tracer, basis_into, basis_from;
	double**	flux, qf[2];
	int			nBasis	= phi->basis[0]->nFuncs;

	flux = new double*[grid->nCells];
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		flux[cell_i] = new double[nBasis];
		for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
			flux[cell_i][basis_i] = 0.0;
		}
	}

	for( edge_i = 0; edge_i < grid->nEdges; edge_i++ ) {
        prePoly = CreatePreImage( edge_i, grid, preGrid, &into, &from, pinds );
		if( prePoly == NULL ) {
			continue;
		}

        for( cell_i = 0; cell_i < 6; cell_i++ ) {
			from = pinds[cell_i];
			incPoly = grid->cells[from];
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

	/* add rhs contributions from previous cell */
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		for( tri_i = 0; tri_i < grid->cells[cell_i]->n; tri_i++ ) {
			tri = grid->cells[cell_i]->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
				weight = tri->wi[quad_i]*tri->Area()/grid->cells[cell_i]->Area();
				tracer = phi->EvalAtCoord( tri->qi[quad_i] );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					basis_into = phi->basis[cell_i]->EvalIJ( qf, basis_i );
					flux[cell_i][basis_i] += weight*tracer*basis_into;
				}
			}
		}

		/* update the cell coefficients */
		AXEB( betaInv_ij[cell_i], flux[cell_i], phiTemp->basis[cell_i]->ci, nBasis );
	}

	if( phiTemp->basis[0]->order == 2 ) {
		Limiter_FirstOrder( phiTemp );
	}
	if( phiTemp->basis[0]->order > 2 ) {
		Limiter_SecondOrder( phiTemp );
	}

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		delete[] flux[cell_i];
	}
	delete[] flux;
}

void CDG::Limiter_FirstOrder( Field* phiTemp ) {
	Grid*	grid	= phi->grid;
	Basis*	basis;
	int		nx, ny;
	int		cell_i, basis_i;
	double	alpha, dxTmp, dyTmp;

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		basis = phiTemp->basis[cell_i];
		nx = cell_i%grid->nx;
		ny = cell_i/grid->nx;

		if( nx > 0 && nx < grid->nx - 1 && ny > 0 && ny < grid->ny - 1 ) {
			alpha = CellLimiter( phiTemp, cell_i );

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
	int		cell_i, basis_i;
	double	alpha_1, alpha_x, alpha_y, alpha_2, dxTmp, dyTmp, dxxTmp, dyyTmp, dxyTmp;

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		basis = phiTemp->basis[cell_i];
		nx = cell_i%grid->nx;
		ny = cell_i/grid->nx;

		if( nx > 0 && nx < grid->nx - 1 && ny > 0 && ny < grid->ny - 1 ) {
			alpha_1 = CellLimiter( phiTemp, cell_i );
			alpha_x = CellSecondLimiter( phiTemp, cell_i, 0 );
			alpha_y = CellSecondLimiter( phiTemp, cell_i, 1 );

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

double CDG::CellLimiter( Field* phiTemp, int cell_i ) {
	Grid*		grid			= phi->grid;
	Polygon* 	poly 			= grid->cells[cell_i];
	Basis*		basis 			= phiTemp->basis[cell_i];
	int			vinds[4];
	int			nx				= cell_i%grid->nx;
	int			ny				= cell_i/grid->nx;
	int			vi, ci, cj;
	double		maxPhiAtVert	= -1.0e+99;
	double		minPhiAtVert	= +1.0e+99;
	double		maxPhiInCell	= -1.0e+99;
	double		minPhiInCell	= +1.0e+99;
	double		phiAtVert, phiInCell, alpha, gamma;

	grid->GetCellVertInds( cell_i, vinds );

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
			phiInCell = phiTemp->basis[cj*grid->nx+ci]->ci[0] - basis->ci[0];
			minPhiInCell = ( phiInCell < minPhiInCell ) ? phiInCell : minPhiInCell;
			maxPhiInCell = ( phiInCell > maxPhiInCell ) ? phiInCell : maxPhiInCell;
		}
	}

	alpha = 1.0;
	if( fabs( minPhiAtVert ) > 1.0e-6 ) {
		alpha = minPhiInCell/minPhiAtVert > 0.0 ? minPhiInCell/minPhiAtVert : 0.0;
	}

	gamma = 1.0;
	if( fabs( maxPhiAtVert ) > 1.0e-6 ) {
		gamma = maxPhiInCell/maxPhiAtVert > 0.0 ? maxPhiInCell/maxPhiAtVert : 0.0;
	}

	alpha = ( alpha < gamma ) ? alpha : gamma;
	alpha = ( alpha < 1.0 )   ? alpha : 1.0;

	return alpha;
}

double CDG::CellSecondLimiter( Field* phiTemp, int cell_i, int dim ) {
	Grid*		grid			= phi->grid;
	Polygon* 	poly 			= grid->cells[cell_i];
	Basis*		basis 			= phiTemp->basis[cell_i];
	int			vinds[4];
	int			nx				= cell_i%grid->nx;
	int			ny				= cell_i/grid->nx;
	int			vi, ci, cj, co;
	double		maxDPhiAtVert	= -1.0e+99;
	double		minDPhiAtVert	= +1.0e+99;
	double		maxDPhiInCell	= -1.0e+99;
	double		minDPhiInCell	= +1.0e+99;
	double		dPhiAtVert, dPhiInCell;
	double 		alpha, gamma;

	grid->GetCellVertInds( cell_i, vinds );

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
			dPhiInCell = phiTemp->basis[co]->EvalDerivFull( grid->cells[co]->origin, dim );
			minDPhiInCell = ( dPhiInCell < minDPhiInCell ) ? dPhiInCell : minDPhiInCell;
			maxDPhiInCell = ( dPhiInCell > maxDPhiInCell ) ? dPhiInCell : maxDPhiInCell;
		}
	}

	alpha = 1.0;
	if( fabs( minDPhiAtVert ) > 1.0e-6 ) {
		alpha = minDPhiInCell/minDPhiAtVert > 0.0 ? minDPhiInCell/minDPhiAtVert : 0.0;
	}

	gamma = 1.0;
	if( fabs( maxDPhiAtVert ) > 1.0e-6 ) {
		gamma = maxDPhiInCell/maxDPhiAtVert > 0.0 ? maxDPhiInCell/maxDPhiAtVert : 0.0;
	}

	alpha = ( alpha < gamma ) ? alpha : gamma;
	alpha = ( alpha < 1.0 )   ? alpha : 1.0;

	return alpha;
}
