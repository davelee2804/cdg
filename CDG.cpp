#include <cstdlib>

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
	int 		i, j, k, l, pi;
	Grid* 		grid 		= phi->grid;
	Cell*		cell;
	Triangle*	tri;
	Basis*		basis;
	int			nBasis		= phi->basis[0]->nFuncs;
	double		fj[nBasis];
	double		weight, *coord;
	double		beta_ij[nBasis*nBasis];

	betaInv_ij = new double*[grid->nCells];

	/* set up the matrix inverse and intial basis coefficients */
	for( pi = 0; pi < grid->nCells; pi++ ) {
		betaInv_ij[pi] = new double[nBasis*nBasis];
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
					fj[j] += weight*basis->EvalIJ( coord, j )*basis->ci[j];
				}
			}
		}
		/* set the initial basis coefficients */
		AXEB( betaInv_ij[pi], fj, basis->ci, nBasis );
	}
}

CDG::~CDG() {
	int i;

	for( i = 0; i < phi->grid->nCells; i++ ) {
		delete[] betaInv_ij[i];
	}
	delete[] betaInv_ij;
}

void CDG::Advect( double dt ) {
	Grid* 	grid	= phi->grid;
	Grid*	preGrid = new Grid( grid->nx, grid->ny, grid->minx, grid->miny, grid->maxx, grid->maxy, grid->quadOrder, grid->basisOrder, true );
	Field*	phiTemp = new Field( grid );

	CalcChars( preGrid, dt );
	preGrid->UpdateEdges();
	preGrid->UpdateCells();
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
	double 		weight, tracer, basis_in, basis_out;
	int			nBasis	= phi->basis[0]->nFuncs;
	double**	flux, origin[2];
	Basis*		basis;

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
				TraceRK2( dt, ADV_BACKWARD, phi->basis[from]->origin, origin );
				basis = new Basis( phi->basis[from]->order, origin );
				//TraceRK2( dt, ADV_BACKWARD, phi->basis[into]->origin, origin );
				//basis = new Basis( phi->basis[into]->order, origin );

				for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
					tri = intPoly->tris[tri_i];
					for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
						weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
						tracer = 0.0;
						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							tracer += phi->basis[from]->ci[basis_i]*basis->EvalIJ( tri->qi[quad_i], basis_i );
							//tracer += phi->basis[into]->ci[basis_i]*basis->EvalIJ( tri->qi[quad_i], basis_i );
						}
						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							basis_in = phi->basis[into]->EvalIJ( tri->qi[quad_i], basis_i );
							basis_out = phi->basis[from]->EvalIJ( tri->qi[quad_i], basis_i );
							flux[into][basis_i] += weight*tracer*basis_in;
							flux[from][basis_i] -= weight*tracer*basis_out;
						}
					}
				}
				delete intPoly;
				delete basis;
			}
        }
        delete prePoly;
	}

	/* add rhs contributions from previous cell */
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		TraceRK2( dt, ADV_BACKWARD, phi->basis[cell_i]->origin, origin );
		basis = new Basis( phi->basis[cell_i]->order, origin );
		for( tri_i = 0; tri_i < grid->cells[cell_i]->n; tri_i++ ) {
			tri = grid->cells[cell_i]->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				weight = tri->wi[quad_i]*tri->Area()/grid->cells[cell_i]->Area();
				tracer = 0.0;
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					tracer += phi->basis[cell_i]->ci[basis_i]*basis->EvalIJ( tri->qi[quad_i], basis_i );
				}
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					basis_in = basis->EvalIJ( tri->qi[quad_i], basis_i );
					flux[cell_i][basis_i] += weight*tracer*basis_in;
				}
			}
		}
		delete basis;

		/* update the cell coefficients */
		AXEB( betaInv_ij[cell_i], flux[cell_i], phiTemp->basis[cell_i]->ci, nBasis );
	}

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		delete[] flux[cell_i];
	}
	delete[] flux;
}

#if 0
void CDG::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 		cell_i, edge_i, basis_i, tri_i, quad_i;
	Grid*		grid	= phi->grid;
	Polygon 	*prePoly, *intPoly;
	Cell		*incPoly;
	Triangle*	tri;
	int			pinds[6], into, from;
	double 		weight, qf[2], tracer, basis_in, basis_out;
	int			order	= grid->basisOrder;
	int			nBasis	= order*order;
	int			nBasis2	= nBasis*nBasis;
	double**	flux;
	double		F_kj[nBasis];
	double		F_kpi[nBasis];
	double		Pij[nBasis2];

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

		/*for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
			F_kj[basis_i] = 0.0;
			F_kpi[basis_i] = 0.0;
		}*/

        for( cell_i = 0; cell_i < 6; cell_i++ ) {
			/*if( pinds[cell_i] == into ) {
				continue;
			}*/

			incPoly = grid->cells[pinds[cell_i]];
			//intPoly = prePoly->Intersection( incPoly );
			intPoly = Intersection( prePoly, incPoly );

			if( intPoly ) {
#ifdef CDG_TEST
				if( pinds[cell_i] == into ) {
					cerr << "ERROR: swept region intersection with inward fluxing cell, area fraction: " << intPoly->Area()/grid->dx/grid->dy << endl;
					abort();
					//continue;
				}
#endif

				/*for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					F_kj[basis_i] = 0.0;
					F_kpi[basis_i] = 0.0;
				}*/

				for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
					tri = intPoly->tris[tri_i];
					for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
						TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );

#ifdef CDG_TEST
						if( qf[0] < grid->minx || qf[0] > grid->maxx || qf[1] < grid->miny || qf[1] > grid->maxy ) {
							cerr << "ERROR: quadrature point outside domain..." << endl;
							incPoly->Print();
							cout << "[" << qf[0] << ", " << qf[1] << "]" << endl; 
							continue;
						}

						if( !intPoly->IsInside( tri->qi[quad_i] ) ) {
							cerr << "ERROR: quadrature point outside intersecting polygon...";
							if( grid->FindCell( tri->qi[quad_i] ) != incPoly ) {
								cerr << " or the incident polygon...";
								abort();
							}
							cerr << endl;
							intPoly->Print();
							cout << "[" << intPoly->verts[0][0] << ", " << intPoly->verts[0][1] << "]" << endl; 
							cout << "[" << tri->qi[quad_i][0] << ", " << tri->qi[quad_i][1] << "]" << endl; 
							continue;
						}
#endif

						weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
						tracer = phi->basis[pinds[cell_i]]->EvalFull( tri->qi[quad_i] );

						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							basis_in = phi->basis[into]->EvalIJ( qf, basis_i );
							basis_out = phi->basis[pinds[cell_i]]->EvalIJ( tri->qi[quad_i], basis_i );

							flux[into][basis_i] += weight*tracer*basis_in;
							flux[pinds[cell_i]][basis_i] -= weight*tracer*basis_out;

							//F_kj[basis_i] -= weight*tracer*basis_out;
						}
					}
				}
				delete intPoly;

				/*BasisProjection( into, pinds[cell_i], Pij );
				AXEB( Pij, F_kj, F_kpi, nBasis );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					flux[into][basis_i] -= F_kpi[basis_i];
					flux[pinds[cell_i]][basis_i] += F_kj[basis_i];
				}*/
			}
        }
        delete prePoly;

		/*BasisProjection( into, from, Pij );
		AXEB( Pij, F_kj, F_kpi, nBasis );
		for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
			flux[into][basis_i] -= F_kpi[basis_i];
			flux[from][basis_i] += F_kj[basis_i];
		}*/
	}
		
	/* add rhs contributions from previous cell */
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		for( tri_i = 0; tri_i < grid->cells[cell_i]->n; tri_i++ ) {
			tri = grid->cells[cell_i]->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );

				weight = tri->wi[quad_i]*tri->Area()/grid->cells[cell_i]->Area();
				//tracer = phi->basis[cell_i]->EvalFull( tri->qi[quad_i] );
				tracer = phi->basis[cell_i]->EvalFull( qf );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					//basis_in = phi->basis[cell_i]->EvalIJ( tri->qi[quad_i], basis_i );
					basis_in = phi->basis[cell_i]->EvalIJ( qf, basis_i );
					flux[cell_i][basis_i] += weight*tracer*basis_in;
				}
			}
		}

		/* update the cell coefficients */
		AXEB( betaInv_ij[cell_i], flux[cell_i], phiTemp->basis[cell_i]->ci, nBasis );
	}

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		delete[] flux[cell_i];
	}
	delete[] flux;
}
#endif
