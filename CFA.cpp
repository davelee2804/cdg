#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
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

#define MAX_POLY_SIZE 12

CFA::CFA( Field* _phi, Field* _velx, Field* _vely, Func* _fu, Func* _fv ) {
	int i;

	phi  = _phi;
	velx = _velx;
	vely = _vely;
	fu   = _fu;
	fv   = _fv;

	pts = new double*[MAX_POLY_SIZE];
	for( i = 0; i < MAX_POLY_SIZE; i++ ) {
		pts[i] = new double[2];
	}
}

CFA::~CFA() {
	int i;

	for( i = 0; i < MAX_POLY_SIZE; i++ ) {
		delete[] pts[i];
	}
	delete[] pts;
}

void CFA::Advect( double dt ) {
	int     i, j;
	Grid*   grid	= phi->grid;
	Grid*   preGrid = new Grid( grid->nx, grid->ny, grid->minx, grid->miny, grid->maxx, grid->maxy, grid->quadOrder, grid->basisOrder, grid->internal );
	Field*  phiTemp = new Field( grid );

	CalcChars( preGrid, dt );
	preGrid->UpdateEdges();
	preGrid->UpdatePolys();
	preGrid->UpdateTris();
	CalcFluxes( preGrid, phiTemp, dt );
	for( i = 0; i < grid->nPolys; i++ ) {
		for( j = 0; j < phi->basis[i]->nFuncs; j++ ) {
			phi->basis[i]->ci[j] += phiTemp->basis[i]->ci[j];
		}
	}

	delete phiTemp;
	delete preGrid;
}

void CFA::TraceEuler( double dt, int dir, double* xi, double* xf ) {
	double xorig[2], vorig[2];

	xorig[0] = xi[0];
	xorig[1] = xi[1];
	if( velx && vely ) {
		velx->LinearInterp( xorig, vorig + 0 );
		vely->LinearInterp( xorig, vorig + 1 );
	}
	else {
		vorig[0] = fu( xorig );
		vorig[1] = fv( xorig );
	}

	/* assume cfl < 1.0 */
	xf[0] = xorig[0] + dir*dt*vorig[0];
	xf[1] = xorig[1] + dir*dt*vorig[1];
	CheckBounds( xf );
}

void CFA::TraceRK2( double dt, int dir, double* xi, double* xf ) {
	double xorig[2], vorig[2], xhalf[2], vhalf[2];

	xorig[0] = xi[0];
	xorig[1] = xi[1];
	if( velx && vely ) {
		velx->LinearInterp( xorig, vorig + 0 );
		vely->LinearInterp( xorig, vorig + 1 );
	}
	else {
		vorig[0] = fu( xorig );
		vorig[1] = fv( xorig );
	}

	/* assume cfl < 1.0 */
	xhalf[0] = xorig[0] + dir*dt*vorig[0];
	xhalf[1] = xorig[1] + dir*dt*vorig[1];
	CheckBounds( xhalf );
	if( velx && vely ) {
		velx->LinearInterp( xhalf, vhalf + 0 );
		vely->LinearInterp( xhalf, vhalf + 1 );
	}
	else {
		vhalf[0] = fu( xhalf );
		vhalf[1] = fv( xhalf );
	}
	xf[0] = xorig[0] + dir*0.5*dt*( vhalf[0] + vorig[0] );
	xf[1] = xorig[1] + dir*0.5*dt*( vhalf[1] + vorig[1] );
	CheckBounds( xf );
}

void CFA::TraceRK4( double dt, int dir, double* xi, double* xf ) {
	double	k1[2], k2[2], k3[2], k4[2], xp[2];

	if( fu == NULL || fv == NULL ) {
		cout << "ERROR: velocity field not supplied as analytic function, so RK4 is not available." << endl;
		abort();
	}

	k1[0] = fu( xi );
	k1[1] = fv( xi );

	xp[0] = xi[0] + dir*0.5*dt*k1[0];
	xp[1] = xi[1] + dir*0.5*dt*k1[1];
	k2[0] = fu( xp );
	k2[1] = fv( xp );

	xp[0] = xi[0] + dir*0.5*dt*k2[0];
	xp[1] = xi[1] + dir*0.5*dt*k2[1];
	k3[0] = fu( xp );
	k3[1] = fv( xp );

	xp[0] = xi[0] + dir*dt*k3[0];
	xp[1] = xi[1] + dir*dt*k3[1];
	k4[0] = fu( xp );
	k4[1] = fv( xp );

	xf[0] = xi[0] + dir*(dt/6.0)*( k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0] );
	xf[1] = xi[1] + dir*(dt/6.0)*( k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1] );

	CheckBounds( xf );
}

void CFA::CalcChars( Grid* preGrid, double dt ) {
	int i, xi, yj;
	int	nx = preGrid->nx;
	int	ny = preGrid->ny;

	#pragma omp parallel private( i, xi, yj )
	{
		#pragma omp for
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
}

double CFA::GetNorm( double* a, double* b, double* c ) {
	double ab[2], ac[2];

	ab[0] = b[0] - a[0];
	ab[1] = b[1] - a[1];
	ac[0] = c[0] - a[0];
	ac[1] = c[1] - a[1];

	return ab[0]*ac[1] - ab[1]*ac[0];
}

void CFA::CheckBounds( double* pt ) {
	if( pt[0] < phi->grid->minx ) pt[0] = phi->grid->minx + 1.0e-8;
	if( pt[0] > phi->grid->maxx ) pt[0] = phi->grid->maxx - 1.0e-8;
	if( pt[1] < phi->grid->miny ) pt[1] = phi->grid->miny + 1.0e-8;
	if( pt[1] > phi->grid->maxy ) pt[1] = phi->grid->maxy - 1.0e-8;
}

Polygon* CFA::CreatePreImage( int ei, Grid* grid, Grid* preGrid, int* into, int* from, int* pinds ) {
	Edge*    e1 = grid->edges[ei];
	Edge*    e2 = preGrid->edges[ei];
	int      left, right, norm;
	Polygon* poly;

	norm = ei/((grid->nx+1)*grid->ny);

	/* ignore boundaries and edges incident on boundaries for now */
	if( !grid->GetEdgePolyInds( ei, pinds ) ) {
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
	//*into = ( GetNorm( grid->edges[ei]->v1, grid->edges[ei]->v2, preGrid->edges[ei]->v2 ) > 0.0 ) ? right : left;
	/* TODO: assuming solid body rotation in normal face evaluation here!!!! */
	int xi, yj;
	grid->EdgeIndexToCoord( ei, &norm, &xi, &yj );
	if( norm == 0 ) {
		*into = ( yj < grid->ny/2 ) ? right : left;
	}
	else {
		*into = ( xi < grid->nx/2 ) ? right : left;
	}
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
	int         ei, pi, pinds[6], into, from;
	Grid*       grid    = phi->grid;
	Polygon     *prePoly, *intPoly, *incPoly;
	double      weight;

	for( ei = 0; ei < grid->nEdges; ei++ ) {
		prePoly = CreatePreImage( ei, grid, preGrid, &into, &from, pinds );
		if( prePoly == NULL ) {
			continue;
		}

		for( pi = 0; pi < 6; pi++ ) {
			from = pinds[pi];
			incPoly = grid->polys[from];
			//intPoly = prePoly->Intersection( incPoly );
			intPoly = Intersection( prePoly, incPoly );
			if( intPoly ) {
#ifdef CFA_TEST
				if( pinds[pi] == into && intPoly->Area()/grid->dx/grid->dy > 1.0e-8 ) {
					cerr << "ERROR: swept region intersection with inward fluxing poly, area fraction: " << intPoly->Area()/grid->dx/grid->dy << endl;
					//abort();
					//continue;
				}
#endif
				weight = intPoly->Area()/incPoly->Area();
				phiTemp->basis[into]->ci[0] += weight*phi->basis[from]->ci[0];
				phiTemp->basis[from]->ci[0] -= weight*phi->basis[from]->ci[0];
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

	/* add points in reverse order for consistency with clockwise polygon convection */
	for( i = 0; i < n; i++ ) {
		pts[i][0] = intersect[n-1-i].x;
		pts[i][1] = intersect[n-1-i].y;
	}

	poly = new Polygon( pts, n, poly1->tris[0]->order );

	return poly;
}
