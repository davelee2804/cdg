#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"
#include "CFA.h"

#include "./core/core"
#include <vector>

//#define CFA_TEST 1

#ifdef CFA_TEST
#include <iostream>
#include <cstdlib>
using namespace std;
#endif

#define ADV_FORWARD (+1)
#define ADV_BACKWARD (-1)

CFA::CFA( Field* _phi, Field* _velx, Field* _vely ) {
	phi  = _phi;
	velx = _velx;
	vely = _vely;

	pts = new double*[4];
	pts[0] = new double[2];
	pts[1] = new double[2];
	pts[2] = new double[2];
	pts[3] = new double[2];
}

CFA::~CFA() {
	delete[] pts[0];
	delete[] pts[1];
	delete[] pts[2];
	delete[] pts[3];
	delete[] pts;
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

Polygon* CFA::CreatePreImage( int ei, Grid* grid, Grid* preGrid, int* into, int* from, int* pinds ) {
	Edge* e1 = grid->edges[ei];
	Edge* e2 = preGrid->edges[ei];
	int left, right, norm;
	Polygon* poly;

	norm = ei/((grid->nx+1)*grid->ny);

	/* ignore boundaries and edges incident on boundaries for now */
	if( !grid->GetEdgeCellInds( ei, pinds ) ) {
		return NULL;
	}

	if( norm == 0 ) {
		left = pinds[2];
		right = pinds[3];

		/* verts must be clockwise */
		pts[0][0] = e1->v2[0];
		pts[0][1] = e1->v2[1];
		pts[1][0] = e1->v1[0];
		pts[1][1] = e1->v1[1];
		pts[2][0] = e2->v1[0];
		pts[2][1] = e2->v1[1];
		pts[3][0] = e2->v2[0];
		pts[3][1] = e2->v2[1];
	}
	else {
		left = pinds[4];
		right = pinds[1];

		/* verts must be clockwise */
		pts[0][0] = e1->v1[0];
		pts[0][1] = e1->v1[1];
		pts[1][0] = e1->v2[0];
		pts[1][1] = e1->v2[1];
		pts[2][0] = e2->v2[0];
		pts[2][1] = e2->v2[1];
		pts[3][0] = e2->v1[0];
		pts[3][1] = e2->v1[1];
	}

	poly = new Polygon( pts, 4, preGrid->quadOrder );

	/* if the cross product of the edge and the vector made by the lower point of the original edge and the 
	   upper point of the final edge is > 0, then the flux is rightwards across the edge */
	*into = ( GetNorm( grid->edges[ei]->v1, grid->edges[ei]->v2, preGrid->edges[ei]->v2 ) > 0.0 ) ? right : left;
	*from = ( *into == left ) ? right : left;

#ifdef CFA_TEST
	Edge* te1 = new Edge( pts[0], pts[3] );
	Edge* te2 = new Edge( pts[1], pts[2] );
	double pt[2];
	if( te1->Intersection( te2, pt ) ) {
		cerr << "ERROR: swept region is a bowtie, edge id:" << ei << endl; 
		abort();
	}
	delete te1;
	delete te2;
#endif

	return poly;
}

void CFA::CalcFluxes( Grid* preGrid, Field* phiTemp, double dt ) {
	int 	ei, pi, pinds[6], into, from;
	Grid*	grid	= phi->grid;
	Polygon	*prePoly, *intPoly;
	Cell*	incPoly;
	double 	weight;

	for( ei = 0; ei < grid->nEdges; ei++ ) {
		prePoly = CreatePreImage( ei, grid, preGrid, &into, &from, pinds );
		if( prePoly == NULL ) {
			continue;
		}

		for( pi = 0; pi < 6; pi++ ) {
			incPoly = grid->cells[pinds[pi]];
			//intPoly = prePoly->Intersection( incPoly );
			intPoly = Intersection( prePoly, incPoly );
			if( intPoly ) {
#ifdef CFA_TEST
				if( pinds[pi] == into && intPoly->Area()/grid->dx/grid->dy > 1.0e-8 ) {
					cerr << "ERROR: swept region intersection with inward fluxing cell, area fraction: " << intPoly->Area()/grid->dx/grid->dy << endl;
					//abort();
					//continue;
				}
#endif
				weight = intPoly->Area()/incPoly->Area();
				phiTemp->basis[into]->ci[0] += weight*phi->basis[pinds[pi]]->ci[0];
				phiTemp->basis[pinds[pi]]->ci[0] -= weight*phi->basis[pinds[pi]]->ci[0];
				delete intPoly;
			}
		}
		delete prePoly;
	}
}

Polygon* CFA::Intersection( Polygon* poly1, Polygon* poly2 ) {
	int i, n;

	typedef core::point<double> Point;
	typedef std::vector<Point> Vec_Point;

	Vec_Point p1(poly1->n), p2(poly2->n), intersect;
	double** pts;
	Polygon* poly;

	for( i = 0; i < poly1->n; i++ ) {
		p1[i].x = poly1->verts[i][0];
		p1[i].y = poly1->verts[i][1];
	}
	for( i = 0; i < poly2->n; i++ ) {
		p2[i].x = poly2->verts[i][0];
		p2[i].y = poly2->verts[i][1];
	}

	core::polygon_intersection( p1, p2, intersect );

	n = intersect.size();
	if( !n ) {
		return NULL;
	}

	pts = new double*[n];
	for( i = 0; i < n; i++ ) {
		pts[i] = new double[2];
		pts[i][0] = intersect[n-1-i].x;
		pts[i][1] = intersect[n-1-i].y;
	}

	poly = new Polygon( pts, n, poly1->tris[0]->order );

	for( i = 0; i < n; i++ ) {
		delete[] pts[i];
	}
	delete[] pts;

	return poly;
}
