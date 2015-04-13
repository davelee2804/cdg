#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"
#include "CFA.h"

#define CFA_TEST 1

#ifdef CFA_TEST
#include <iostream>
#include <cstdlib>
using namespace std;
#endif

#define ADV_FORWARD +1
#define ADV_BACKWARD -1

CFA::CFA( Field* _phi, Field* _velx, Field* _vely ) {
	phi  = _phi;
	velx = _velx;
	vely = _vely;
}

CFA::~CFA() {
}

void CFA::Advect( double dt ) {
	int		i, j;
	Grid* 	grid	= phi->grid;
	Grid*	preGrid = new Grid( grid->nx, grid->ny, grid->minx, grid->miny, grid->maxx, grid->maxy, grid->quadOrder, grid->basisOrder, grid->internal );
	Field*	phiTemp = new Field( grid );

	CalcChars( preGrid, dt );
	preGrid->UpdateEdges();
	preGrid->UpdateCells();
	CalcFluxes( preGrid, phiTemp, dt );
	for( i = 0; i < grid->nCells; i++ ) {
		for( j = 0; j < grid->cells[0]->nc; j++ ) {
			phi->basis[i]->ci[j] += phiTemp->basis[i]->ci[j];
		}
	}

	delete phiTemp;
	delete preGrid;
}

void CFA::TraceRK2( double dt, int dir, double* xi, double* xf ) {
	double xorig[2], vorig[2], xhalf[2], vhalf[2];

	xorig[0] = xi[0];
	xorig[1] = xi[1];
	velx->LinearInterp( xorig, vorig + 0 );
	vely->LinearInterp( xorig, vorig + 1 );

	/* assume cfl < 1.0 */
	xhalf[0] = xorig[0] + dir*dt*vorig[0];
	xhalf[1] = xorig[1] + dir*dt*vorig[1];
	velx->LinearInterp( xhalf, vhalf + 0 );
	vely->LinearInterp( xhalf, vhalf + 1 );
	xf[0] = xorig[0] + dir*0.5*dt*( vhalf[0] + vorig[0] );
	xf[1] = xorig[1] + dir*0.5*dt*( vhalf[1] + vorig[1] );
}

void CFA::CalcChars( Grid* preGrid, double dt ) {
	int i, xi, yj;
	int	nx = preGrid->nx;
	int	ny = preGrid->ny;

	/* calculate the grid pre-image */
	for( i = 0; i < phi->grid->nVerts; i++ ) {
		xi = i%(nx+1);
		yj = i/(nx+1);

		/* ignore boundary verts for now */
		if( xi == 0 || xi == nx || yj == 0 || yj == ny ) {
			continue;
		}

		TraceRK2( dt, ADV_BACKWARD, phi->grid->verts[i], preGrid->verts[i] );
	}
}

double CFA::GetNorm( double* a, double* b, double* c ) {
	double ab[2], ac[2];

	ab[0] = b[0] - a[0];
	ab[1] = b[1] - a[1];
	ac[0] = c[0] - a[0];
	ac[1] = c[1] - a[1];

	return ab[0]*ac[1] - ab[1]*ac[0];
}

void CFA::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 	ei, norm, xi, yj, pi;
	Grid*	grid	= phi->grid;
	double 	**pts;
	Polygon	*prePoly, *intPoly;
	Cell*	incPoly;
	int		pinds[6], left, right, into;
	double 	weight;

	pts = new double*[4];
	pts[0] = new double[2];
	pts[1] = new double[2];
	pts[2] = new double[2];
	pts[3] = new double[2];

	for( ei = 0; ei < grid->nEdges; ei++ ) {
		grid->EdgeIndexToCoord( ei, &norm, &xi, &yj );

		/* ignore boundaries and  edges incident on boundaries for now */
		if( !grid->GetEdgeCellInds( ei, pinds ) ) {
			continue;
		}

		if( norm == 0 ) {
			left = pinds[2];
			right = pinds[3];

			/* verts must be clockwise */
			pts[0][0] = grid->edges[ei]->v2[0];
			pts[0][1] = grid->edges[ei]->v2[1];
			pts[1][0] = grid->edges[ei]->v1[0];
			pts[1][1] = grid->edges[ei]->v1[1];
			pts[2][0] = preGrid->edges[ei]->v1[0];
			pts[2][1] = preGrid->edges[ei]->v1[1];
			pts[3][0] = preGrid->edges[ei]->v2[0];
			pts[3][1] = preGrid->edges[ei]->v2[1];
		}
		else {
			left = pinds[4];
			right = pinds[1];

			/* verts must be clockwise */
			pts[0][0] = grid->edges[ei]->v1[0];
			pts[0][1] = grid->edges[ei]->v1[1];
			pts[1][0] = grid->edges[ei]->v2[0];
			pts[1][1] = grid->edges[ei]->v2[1];
			pts[2][0] = preGrid->edges[ei]->v2[0];
			pts[2][1] = preGrid->edges[ei]->v2[1];
			pts[3][0] = preGrid->edges[ei]->v1[0];
			pts[3][1] = preGrid->edges[ei]->v1[1];
		}

#ifdef CFA_TEST
		Edge* e1 = new Edge( pts[0], pts[3] );
		Edge* e2 = new Edge( pts[1], pts[2] );
		double pt[2];
		if( e1->Intersection( e2, pt ) ) {
			cerr << "ERROR: swept region is a bowtie, edge id:" << ei << endl; 
			abort();
		}
		delete e1;
		delete e2;
#endif

		prePoly = new Polygon( pts, 4, preGrid->quadOrder );

		/* if the cross product of the edge and the vector made by the lower point of the original edge and the 
		   upper point of the final edge is > 0, then the flux is rightwards across the edge */
		into = ( GetNorm( grid->edges[ei]->v1, grid->edges[ei]->v2, preGrid->edges[ei]->v2 ) > 0.0 ) ? right : left;
		for( pi = 0; pi < 6; pi++ ) {
			if( pinds[pi] == into ) {
				continue;
			}

			incPoly = grid->cells[pinds[pi]];
			intPoly = prePoly->Intersection( incPoly );
			if( intPoly ) {
				weight = intPoly->Area()/incPoly->Area();
				phiTemp->basis[into]->ci[0] += weight*phi->basis[pinds[pi]]->ci[0];
				phiTemp->basis[pinds[pi]]->ci[0] -= weight*phi->basis[pinds[pi]]->ci[0];
				delete intPoly;
			}
		}
		delete prePoly;
	}

	delete[] pts[0];
	delete[] pts[1];
	delete[] pts[2];
	delete[] pts[3];
	delete[] pts;
}
