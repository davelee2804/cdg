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

CDG::CDG( Field* _phi, Field* _velx, Field* _vely, Func* _fu, Func* _fv ) : CFA( _phi, _velx, _vely, _fu, _fv ) {
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
	int        i, j, k, l, pi;
	Grid*      grid   = phi->grid;
	Polygon*   poly;
	Triangle*  tri;
	Basis*     basis;
	int        nBasis = phi->basis[0]->nFuncs;
	double     fj[nBasis];
	double     weight, *coord;
	double     beta_ij[nBasis*nBasis];
	double     volErr;

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
			for( l = 0; l < tri->nq; l++ ) {
				for( j = 0; j < nBasis; j++ ) {
					for( i = 0; i < nBasis; i++ ) {
						weight = tri->wi[l]*tri->area;
						coord = tri->qi[l];
						beta_ij[j*nBasis+i] += weight*basis->EvalIJ( coord, i )*basis->EvalIJ( coord, j );
					}
				}
			}
		}
		MatInv( beta_ij, betaInv_ij[pi], nBasis );

		/* if we're not supplied a function, then we're not initializing the basis coefficients */
		if( !func ) {
			continue;
		}

		for( j = 0; j < nBasis; j++ ) {
			fj[j] = 0.0;
			for( k = 0; k < poly->n; k++ ) {
				tri = poly->tris[k];
				for( l = 0; l < tri->nq; l++ ) {
					weight = tri->wi[l]*tri->area;
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
	int     i;
	Grid*   grid    = phi->grid;
	Grid*   preGrid = new Grid( grid->nx, grid->ny, grid->minx, grid->miny, grid->maxx, grid->maxy, grid->quadOrder, grid->basisOrder, true );
	Field*  phiTemp = new Field( grid );

	for( i = 0; i < grid->nVerts; i++ ) {
		preGrid->verts[i][0] = grid->verts[i][0];
		preGrid->verts[i][1] = grid->verts[i][1];
	}

	CalcChars( preGrid, dt );
	preGrid->UpdateEdges();
	preGrid->UpdatePolys();
	preGrid->UpdateTris();
	CalcFluxes( preGrid, phiTemp, dt );
	phi->Copy( phiTemp );

	delete phiTemp;
	delete preGrid;
}

void CDG::BasisProjection( int kp, int k, double* Pij ) {
	int         tri_i, quad_i, basis_m, basis_j;
	Grid*       grid        = phi->grid;
	Polygon*    poly        = grid->polys[kp];
	Basis*      basis_k     = phi->basis[k];
	Basis*      basis_kp    = phi->basis[kp];
	Triangle*   tri;
	int         nBasis      = phi->basis[0]->nFuncs;
	int         nBasis2     = nBasis*nBasis;
	double      beta_mj[nBasis2];
	double      weight, *coord;

	for( basis_j = 0; basis_j < nBasis2; basis_j++ ) {
		beta_mj[basis_j] = 0.0;
	}

	for( tri_i = 0; tri_i < poly->n; tri_i++ ) {
		tri = poly->tris[tri_i];
		for( quad_i = 0; quad_i < tri->nq; quad_i++ ) {
			weight = tri->wi[quad_i]*tri->area;
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
	int         poly_i, edge_i, basis_i, tri_i, quad_i;
	Grid*       grid    = phi->grid;
	Polygon     *prePoly, *intPoly, *incPoly;
	Triangle*   tri;
	int         pinds[6], into, from;
	double      weight, tracer, basis_into, basis_from;
	double**    flux, qf[2];
	int         nBasis  = phi->basis[0]->nFuncs;

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
			incPoly = grid->polys[pinds[poly_i]];
			intPoly = Intersection( prePoly, incPoly );
			if( intPoly == NULL ) {
				continue;
			}

			for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
				tri = intPoly->tris[tri_i];
				for( quad_i = 0; quad_i < tri->nq; quad_i++ ) {
					TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
					weight = tri->wi[quad_i]*tri->area;
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
		delete prePoly;
	}

	#pragma omp parallel private( poly_i, tri_i, tri, quad_i, qf, weight, tracer, basis_into )
	{
		#pragma omp for
		/* add rhs contributions from previous poly */
		for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
			for( tri_i = 0; tri_i < grid->polys[poly_i]->n; tri_i++ ) {
				tri = grid->polys[poly_i]->tris[tri_i];
				for( quad_i = 0; quad_i < tri->nq; quad_i++ ) {
					TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
					weight = tri->wi[quad_i]*tri->area;
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
	}

	for( poly_i = 0; poly_i < grid->nPolys; poly_i++ ) {
		delete[] flux[poly_i];
	}
	delete[] flux;
}
