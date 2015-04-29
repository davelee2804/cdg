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

	phiMax = -1.0e+10;
	phiMin = +1.0e+10;
	for( pi = 0; pi < grid->nCells; pi++ ) {
		for( j = 0; j < nBasis; j++ ) {
			if( phi->basis[pi]->ci[j] > phiMax ) {
				phiMax = phi->basis[pi]->ci[j];
			}
			if( phi->basis[pi]->ci[j] < phiMin ) {
				phiMin = phi->basis[pi]->ci[j];
			}
		}
	}

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

		if( !basis->TestMean( &volErr ) ) {
			cout << "ERROR: basis function mean not equal to first component..." << volErr << endl;
		}
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
	double**	flux, origin[2];
	Basis*		basis;
	int			nBasis	= phi->basis[0]->nFuncs;
	int			nBasis2	= nBasis*nBasis;
	double		F_kj[nBasis], F_kpi[nBasis], P_ij[nBasis2];

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
				basis = new Basis( preGrid->cells[from], phi->basis[from]->order, origin, grid->dx, grid->dy );

				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					F_kj[basis_i] = 0.0;
					F_kpi[basis_i] = 0.0;
				}

				for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
					tri = intPoly->tris[tri_i];
					for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
						weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
						tracer = basis->EvalWithCoeffs( tri->qi[quad_i], phi->basis[from]->ci );
						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							basis_into = phi->basis[into]->EvalIJ( tri->qi[quad_i], basis_i );
							basis_from = phi->basis[from]->EvalIJ( tri->qi[quad_i], basis_i );
							//flux[into][basis_i] += weight*tracer*basis_into;
							//flux[from][basis_i] -= weight*tracer*basis_from;

							F_kj[basis_i] -= weight*tracer*basis_from;
						}
					}
				}
				delete intPoly;
				delete basis;

				BasisProjection( into, from, P_ij );
				AXEB( P_ij, F_kj, F_kpi, nBasis );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					flux[into][basis_i] -= F_kpi[basis_i];
					flux[from][basis_i] += F_kj[basis_i];
				}
			}
        }
        delete prePoly;
	}

	/* add rhs contributions from previous cell */
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		TraceRK2( dt, ADV_BACKWARD, phi->basis[cell_i]->origin, origin );
		basis = new Basis( preGrid->cells[cell_i], phi->basis[cell_i]->order, origin, grid->dx, grid->dy );
		for( tri_i = 0; tri_i < grid->cells[cell_i]->n; tri_i++ ) {
			tri = grid->cells[cell_i]->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				weight = tri->wi[quad_i]*tri->Area()/grid->cells[cell_i]->Area();
				tracer = basis->EvalWithCoeffs( tri->qi[quad_i], phi->basis[cell_i]->ci );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					basis_into = basis->EvalIJ( tri->qi[quad_i], basis_i );
					flux[cell_i][basis_i] += weight*tracer*basis_into;
				}
			}
		}

		/* update the cell coefficients */
		AXEB( betaInv_ij[cell_i], flux[cell_i], phiTemp->basis[cell_i]->ci, nBasis );
		/* apply the limiter */
		delete basis;
	}

	Limiter( phiTemp );

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		delete[] flux[cell_i];
	}
	delete[] flux;
}

void CDG::Limiter( Field* phiTemp ) {
	Grid*	grid	= phi->grid;
	Cell*	cell;
	Basis*	basis;
	int		nx, ny;
	int		vinds[4];
	int		cell_i, vert_i, ci, cj, basis_i;
	double	maxPhiAtVert, minPhiAtVert, phiAtVert, maxPhiInCell, minPhiInCell, phiInCell;
	double	alpha, gamma, dxTmp, dyTmp;

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		cell = grid->cells[cell_i];
		basis = phiTemp->basis[cell_i];
		nx = cell_i%grid->nx;
		ny = cell_i/grid->nx;

		if( nx > 0 && nx < grid->nx - 1 && ny > 0 && ny < grid->ny - 1 ) {
			grid->GetCellVertInds( cell_i, vinds );
			maxPhiAtVert = maxPhiInCell = -1.0e+99;
			minPhiAtVert = minPhiInCell = +1.0e+99;
			for( vert_i = 0; vert_i < cell->n; vert_i++ ) {
				phiAtVert = basis->EvalFull( cell->verts[vert_i] ) - basis->ci[0];
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

			if( alpha > 1.0 - 1.0e-4 ) {
				continue;
			}

			dxTmp = basis->ci[1];
			dyTmp = basis->ci[basis->order];

			for( basis_i = 1; basis_i < basis->nFuncs; basis_i++ ) {
				basis->ci[basis_i] = 0.0;
			}
			basis->ci[1] = alpha*dxTmp;
			basis->ci[basis->order] = alpha*dyTmp;
		}
	}
}

#if 0
void CDG::Limiter( Field* phiTemp ) {
	Grid*	grid	= phi->grid;
	Cell*	cell;
	Basis*	basis;
	int 	cell_i, vert_i, cell_j, adjCell_i, basis_j;
	int		nx, ny;
	int		vinds[4], cinds[4], adjCells[3], maxInds[2];
	double	atVerts[4], atVerts2[2];
	bool	minLimit, maxLimit;
	double	Aij[4], cj[2], AijInv[4];
	double	incVal[4][3], min, max, avg;

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		cell = grid->cells[cell_i];
		basis = phiTemp->basis[cell_i];
		avg = basis->ci[0];
		minLimit = maxLimit = false;
		nx = cell_i%grid->nx;
		ny = cell_i/grid->nx;

		if( nx > 0 && nx < grid->nx - 1 && ny > 0 && ny < grid->ny - 1 ) {
			/* for each vertex of the cell get the adjacent cell mean values */
			grid->GetCellVertInds( cell_i, vinds );
			for( vert_i = 0; vert_i < cell->n; vert_i++ ) {
				atVerts[vert_i] = basis->EvalFull( cell->verts[vert_i] );
				grid->GetVertCellInds( vinds[vert_i], cinds );
				adjCell_i = 0;
				for( cell_j = 0; cell_j < 4; cell_j++ ) {
					if( cinds[cell_j] == cell_i ) {
						continue;
					}
					adjCells[adjCell_i++] = cinds[cell_j];
				}

				for( adjCell_i = 0; adjCell_i < 3; adjCell_i++ ) {
					incVal[vert_i][adjCell_i] = phiTemp->basis[adjCells[adjCell_i]]->ci[0];
				}
			}

			/* find the highest minimum and lowest maximum  mean cell values exceeded at the vertices */
			min = -1.0e+99;
			max = +1.0e+99;
			for( vert_i = 0; vert_i < cell->n; vert_i++ ) {
				for( adjCell_i = 0; adjCell_i < 3; adjCell_i++ ) {
					/* highest minimum */
					if( atVerts[vert_i] < avg && atVerts[vert_i] < incVal[vert_i][adjCell_i] ) {
						minLimit = true;
						if( incVal[vert_i][adjCell_i] > min ) {
							min = incVal[vert_i][adjCell_i];
							maxInds[0] = vert_i;
						}
					}
					/* lowest maximum */
					else if( atVerts[vert_i] > avg && atVerts[vert_i] > incVal[vert_i][adjCell_i] ) {
						maxLimit = true;
						if( incVal[vert_i][adjCell_i] < max ) {
							max = incVal[vert_i][adjCell_i];
							maxInds[1] = vert_i;
						}
					}
				}
			}

			if( !minLimit && !maxLimit ) {
				continue;
			}

			if( !minLimit ) {
				maxInds[0] = (maxInds[1]+2)%4;
				min = atVerts[maxInds[0]];
			}
			else if( !maxLimit ) {
				maxInds[1] = (maxInds[0]+2)%4;
				max = atVerts[maxInds[1]];
			}

			atVerts2[0] = min;
			atVerts2[1] = max;
			for( vert_i = 0; vert_i < 2; vert_i++ ) {
				for( basis_j = 0; basis_j < 2; basis_j++ ) {
					Aij[vert_i*2+basis_j] = basis->EvalIJ( cell->verts[maxInds[vert_i]], basis_j + 1 );
				}
			}
			MatInv( Aij, AijInv, 2 );
			AXEB( AijInv, atVerts2, cj, 2 );
			for( basis_j = 1; basis_j < basis->nFuncs; basis_j++ ) {
				basis->ci[basis_j] = 0.0;
			}
			basis->ci[1] = cj[1];
			basis->ci[basis->order] = cj[2];
		}
	}
}
#endif

#if 0
void CDG::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 		cell_i, edge_i, basis_i, tri_i, quad_i;
	Grid*		grid	= phi->grid;
	Polygon 	*prePoly, *intPoly;
	Cell		*incPoly;
	Triangle*	tri;
	int			pinds[6], into, from;
	double 		weight, qf[2], tracer, basis_into, basis_from;
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
							basis_into = phi->basis[into]->EvalIJ( qf, basis_i );
							basis_from = phi->basis[pinds[cell_i]]->EvalIJ( tri->qi[quad_i], basis_i );

							flux[into][basis_i] += weight*tracer*basis_into;
							flux[pinds[cell_i]][basis_i] -= weight*tracer*basis_from;

							//F_kj[basis_i] -= weight*tracer*basis_from;
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
					//basis_into = phi->basis[cell_i]->EvalIJ( tri->qi[quad_i], basis_i );
					basis_into = phi->basis[cell_i]->EvalIJ( qf, basis_i );
					flux[cell_i][basis_i] += weight*tracer*basis_into;
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
