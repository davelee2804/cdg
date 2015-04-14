#include "Edge.h"
#include "Triangle.h"
#include "Basis.h"
#include "Polygon.h"
#include "Cell.h"
#include "Grid.h"
#include "Field.h"
#include "CFA.h"
#include "CDG.h"

#define ADV_FORWARD +1
#define ADV_BACKWARD -1

CDG::CDG( Field* _phi, Field* _velx, Field* _vely ) : CFA( _phi, _velx, _vely ) {
	int 		i, j, k, l, pi;
	Grid* 		grid 		= phi->grid;
	Cell*		cell;
	Triangle*	tri;
	Basis*		basis;
	int			order		= phi->basis[0]->order;
	int			nBasis		= order*order;
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
						//beta_ij[j*nBasis+i] += weight*basis->EvalIJ( cell->coords[i], i, cell->coords[i] )*basis->EvalIJ( cell->coords[j], j, cell->coords[j] );
						beta_ij[j*nBasis+i] += weight*basis->EvalIJ( cell->coords[i], i, coord )*basis->EvalIJ( cell->coords[j], j, coord );
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
					/* basis initially set as the spatial values at the cell coordinates */
					fj[j] += weight*basis->EvalIJ( cell->coords[j], j, cell->coords[j] )*basis->ci[j];
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

void CDG::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 		cell_i, cell_x, cell_y, edge_i, basis_i, tri_i, quad_i, norm;
	Grid*		grid	= phi->grid;
	int			nx 		= grid->nx;
	int			ny 		= grid->ny;
	double 		**pts;
	Polygon 	*prePoly, *intPoly;
	Cell		*incPoly;
	Triangle*	tri;
	int			pinds[6], left, right, into;
	double 		weight, qf[2], tracer, basis, sign;
	int			order	= grid->basisOrder;
	int			nBasis	= order*order;
	double**	flux;
	double**	preCoords;

	pts = new double*[4];
	pts[0] = new double[2];
	pts[1] = new double[2];
	pts[2] = new double[2];
	pts[3] = new double[2];

	flux = new double*[grid->nCells];
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		flux[cell_i] = new double[nBasis];
		for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
			flux[cell_i][basis_i] = 0.0;
		}
	}

	preCoords = new double*[nBasis];
	for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
		preCoords[basis_i] = new double[2];
	}

	for( edge_i = 0; edge_i < grid->nEdges; edge_i++ ) {
        grid->EdgeIndexToCoord( edge_i, &norm, &cell_x, &cell_y );

        /* ignore boundaries and edges incident on boundaries for now */
        if( !grid->GetEdgeCellInds( edge_i, pinds ) ) {
            continue;
        }

        if( norm == 0 ) {
            left = pinds[2];
            right = pinds[3];

            /* verts must be clockwise */
            pts[0][0] = grid->edges[edge_i]->v2[0];
            pts[0][1] = grid->edges[edge_i]->v2[1];
            pts[1][0] = grid->edges[edge_i]->v1[0];
            pts[1][1] = grid->edges[edge_i]->v1[1];
            pts[2][0] = preGrid->edges[edge_i]->v1[0];
            pts[2][1] = preGrid->edges[edge_i]->v1[1];
            pts[3][0] = preGrid->edges[edge_i]->v2[0];
            pts[3][1] = preGrid->edges[edge_i]->v2[1];
        }
        else {
            left = pinds[4];
            right = pinds[1];

            /* verts must be clockwise */
            pts[0][0] = grid->edges[edge_i]->v1[0];
            pts[0][1] = grid->edges[edge_i]->v1[1];
            pts[1][0] = grid->edges[edge_i]->v2[0];
            pts[1][1] = grid->edges[edge_i]->v2[1];
            pts[2][0] = preGrid->edges[edge_i]->v2[0];
            pts[2][1] = preGrid->edges[edge_i]->v2[1];
            pts[3][0] = preGrid->edges[edge_i]->v1[0];
            pts[3][1] = preGrid->edges[edge_i]->v1[1];
        }

        prePoly = new Polygon( pts, 4, preGrid->quadOrder );

        /* if the cross product of the edge and the vector made by the lower point of the original edge and the 
		   upper point of the final edge is > 0, then the flux is rightwards across the edge */
        into = ( GetNorm( grid->edges[edge_i]->v1, grid->edges[edge_i]->v2, preGrid->edges[edge_i]->v2 ) > 0.0 ) ? right : left;
        for( cell_i = 0; cell_i < 6; cell_i++ ) {
			sign = ( pinds[cell_i] == into ) ? +1.0 : -1.0;
			incPoly = grid->cells[pinds[cell_i]];
			intPoly = prePoly->Intersection( incPoly );
			if( intPoly ) {
				for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
					tri = intPoly->tris[tri_i];
					for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
						//TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							TraceRK2( dt, ADV_BACKWARD, grid->cells[pinds[cell_i]]->coords[basis_i], preCoords[basis_i] );
						}

						weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();

						//tracer = phi->basis[pinds[cell_i]]->EvalFull( tri->qi[quad_i], grid->cells[pinds[cell_i]]->coords ); //where to evaluate??
						//tracer = phi->basis[into]->EvalFull( qf, grid->cells[into]->coords ); //where to evaluate??
						//tracer = phi->basis[into]->EvalFull( tri->qi[quad_i], preCoords ); //where to evaluate??
						tracer = phi->basis[pinds[cell_i]]->EvalFull( tri->qi[quad_i], preCoords ); //where to evaluate??

						for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
							//basis = phi->basis[into]->EvalIJ( qf, basis_i, grid->cells[into]->coords[basis_i] );
							basis = phi->basis[into]->EvalIJ( tri->qi[quad_i], basis_i, preCoords[basis_i] );

							flux[pinds[cell_i]][basis_i] += sign*weight*tracer*basis;
						}
					}
				}
				delete intPoly;
			}
        }
        delete prePoly;
	}
		
	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		cell_x = cell_i%nx;
		cell_y = cell_i/nx;

		/* ignore boundaries for now */
		if( cell_x == 0 || cell_x == nx - 1 || cell_y == 0 || cell_y == ny - 1 ) {
			continue;
		}

		/* add rhs contributions from previous cell */
		for( tri_i = 0; tri_i < grid->cells[cell_i]->n; tri_i++ ) {
			tri = grid->cells[cell_i]->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				weight = tri->wi[quad_i]*tri->Area()/grid->cells[cell_i]->Area();
				tracer = phi->basis[cell_i]->EvalFull( tri->qi[quad_i], grid->cells[cell_i]->coords );
				for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
					basis = phi->basis[cell_i]->EvalIJ( tri->qi[quad_i], basis_i, grid->cells[cell_i]->coords[basis_i] );
					flux[cell_i][basis_i] += weight*tracer*basis;
				}
			}
		}

		/* update the cell coefficients */
		AXEB( betaInv_ij[cell_i], flux[cell_i], phiTemp->basis[cell_i]->ci, nBasis );
	}

	delete[] pts[0];
	delete[] pts[1];
	delete[] pts[2];
	delete[] pts[3];
	delete[] pts;

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		delete[] flux[cell_i];
	}
	delete[] flux;

	for( basis_i = 0; basis_i < nBasis; basis_i++ ) {
		delete[] preCoords[basis_i];
	}
	delete[] preCoords;
}

#if 0
void CDG::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 		cell_i, cell_x, cell_y, edge_i, pi, tri_i, quad_i, fj, v1, v2;
	Grid*		grid	= phi->grid;
	int			nx 		= grid->nx;
	int			ny 		= grid->ny;
	double 		**pts;
	Polygon 	*prePoly, *intPoly;
	Cell		*incPoly, *locPoly;
	Triangle*	tri;
	int			pinds[6], einds[4], left, right, into;
	double 		weight, qf[2], tracer, val;
	int			order	= grid->basisOrder;
	int			nBasis	= order*order;
	double		flux[nBasis];

	pts = new double*[4];
	pts[0] = new double[2];
	pts[1] = new double[2];
	pts[2] = new double[2];
	pts[3] = new double[2];

	for( cell_i = 0; cell_i < grid->nCells; cell_i++ ) {
		cell_x = cell_i%nx;
		cell_y = cell_i/nx;

		/* ignore boundaries for now */
		if( cell_x == 0 || cell_x == nx - 1 || cell_y == 0 || cell_y == ny - 1 ) {
			continue;
		}

		grid->GetCellEdgeInds( cell_i, einds );
		locPoly = grid->cells[cell_i];

		for( pi = 0; pi < nBasis; pi++ ) {
			flux[pi] = 0.0;
		}

		for( edge_i = 0; edge_i < 4; edge_i++ ) {
			if( edge_i == 0 ) {
				v1 = (cell_y+1)*(nx+1)+(cell_x+0);
				v2 = (cell_y+0)*(nx+1)+(cell_x+0);
				left  = (cell_y+0)*nx+(cell_x-1);
				right = (cell_y+0)*nx+(cell_x+0);
			}
			else if( edge_i == 1 ) {
				v1 = (cell_y+1)*(nx+1)+(cell_x+0);
				v2 = (cell_y+1)*(nx+1)+(cell_x+1);
				left  = (cell_y+1)*nx+(cell_x+0);
				right = (cell_y+0)*nx+(cell_x+0);
			}
			else if( edge_i == 2 ) {
				v1 = (cell_y+1)*(nx+1)+(cell_x+1);
				v2 = (cell_y+0)*(nx+1)+(cell_x+1);
				left  = (cell_y+0)*nx+(cell_x+0);
				right = (cell_y+0)*nx+(cell_x+1);
			}
			else if( edge_i == 3 ) {
				v1 = (cell_y+0)*(nx+1)+(cell_x+0);
				v2 = (cell_y+0)*(nx+1)+(cell_x+1);
				left  = (cell_y+0)*nx+(cell_x+0);
				right = (cell_y-1)*nx+(cell_x+0);
			}

			/* verts must be clockwise */
			pts[0][0] = grid->verts[v1][0];
			pts[0][1] = grid->verts[v1][1];
			pts[1][0] = grid->verts[v2][0];
			pts[1][1] = grid->verts[v2][1];
			pts[2][0] = preGrid->verts[v2][0];
			pts[2][1] = preGrid->verts[v2][1];
			pts[3][0] = preGrid->verts[v1][0];
			pts[3][1] = preGrid->verts[v1][1];

			prePoly = new Polygon( pts, 4, preGrid->quadOrder );
			into = ( GetNorm( grid->edges[einds[edge_i]]->v1, grid->edges[einds[edge_i]]->v2, preGrid->edges[einds[edge_i]]->v2 ) > 0.0 ) ? right : left;

			/* flux is into cell */
			if( into == cell_i ) {
				grid->GetEdgeCellInds( einds[edge_i], pinds );
				for( pi = 0; pi < 6; pi++ ) {
					incPoly = grid->cells[pinds[pi]];
					intPoly = prePoly->Intersection( incPoly );
					if( intPoly ) {
						for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
							tri = intPoly->tris[tri_i];
							for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
								TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
								weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
								tracer = phi->basis[pinds[pi]]->EvalFull( tri->qi[quad_i], grid->cells[pinds[pi]]->coords );
								for( fj = 0; fj < nBasis; fj++ ) {
									val = phi->basis[cell_i]->EvalIJ( qf, fj, grid->cells[cell_i]->coords[fj] );
									flux[fj] += weight*tracer*val;
								}
							}
						}
						delete intPoly;
					}
				}
			}
			/* flux is out of cell */
			else {
				intPoly = prePoly->Intersection( locPoly );
				if( intPoly ) {
					for( tri_i = 0; tri_i < intPoly->n; tri_i++ ) {
						tri = intPoly->tris[tri_i];
						for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
							TraceRK2( dt, ADV_FORWARD, tri->qi[quad_i], qf );
							into = grid->GetCellIndex( qf );
							weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
							tracer = phi->basis[into]->EvalFull( tri->qi[quad_i], grid->cells[into]->coords );
							for( fj = 0; fj < nBasis; fj++ ) {
								val = phi->basis[into]->EvalIJ( qf, fj, grid->cells[into]->coords[fj] );
								flux[fj] -= weight*tracer*val;
							}
						}
					}
					delete intPoly;
				}
			}
			delete prePoly;
		}

		/* add rhs contributions from previous cell */
		for( tri_i = 0; tri_i < locPoly->n; tri_i++ ) {
			tri = locPoly->tris[tri_i];
			for( quad_i = 0; quad_i < tri->nQuadPts; quad_i++ ) {
				weight = tri->wi[quad_i]*tri->Area()/incPoly->Area();
				tracer = phi->basis[cell_i]->EvalFull( tri->qi[quad_i], grid->cells[cell_i]->coords );
				for( fj = 0; fj < nBasis; fj++ ) {
					val = phi->basis[cell_i]->EvalIJ( tri->qi[quad_i], fj, grid->cells[cell_i]->coords[fj] );
					flux[fj] += weight*tracer*val;
				}
			}
		}

		/* update the cell coefficients */
		AXEB( betaInv_ij[cell_i], flux, phiTemp->basis[cell_i]->ci, nBasis );
	}

	delete[] pts[0];
	delete[] pts[1];
	delete[] pts[2];
	delete[] pts[3];
	delete[] pts;
}
#endif
