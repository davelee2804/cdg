#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"

Polygon::Polygon( double** _verts, int _n, int quadOrder ) {
	int i;

	n = _n;

	verts = new double*[n];
	edges = new Edge*[n];
	tris = new Triangle*[n];

	for( i = 0; i < n; i++ ) {
		verts[i] = new double[2];
		verts[i][0] = _verts[i][0];
		verts[i][1] = _verts[i][1];
	}

	GenOrigin();

	for( i = 0; i < n; i++ ) {
		edges[i] = new Edge( verts[i], verts[(i+1)%n] );
		tris[i] = new Triangle( verts[i], verts[(i+1)%n], origin, quadOrder );
	}
}

Polygon::~Polygon() {
	int i;

	for( i = 0; i < n; i++ ) {
		delete[] verts[i];
		delete edges[i];
		delete tris[i];
	}
	delete[] verts;
	delete[] edges;
	delete[] tris;
}

#define R2_TOL 1.0e-8

int Polygon::AddPoint( double** pts, int np, double* pt ) {
	int i;

	for( i = 0; i < np; i++ ) {
		if( (pts[i][0] - pt[0])*(pts[i][0] - pt[0]) + (pts[i][1] - pt[1])*(pts[i][1] - pt[1]) < R2_TOL ) {
			return np;
		}
	}

	pts[np][0] = pt[0];
	pts[np][1] = pt[1];
	return np + 1;
}

#define MAX_POLY 10

Polygon* Polygon::Intersection( Polygon* poly ) {
	int         vi, vj, iOffset = -1, jOffset = -1, pSize = 0;
	double      **polyVerts;
	double      ijInt[2], jiInt[2];
	Polygon*    intPoly	= NULL;
	bool        allInside;

	polyVerts = new double*[MAX_POLY];
	for( vi = 0; vi < MAX_POLY; vi++ ) {
		polyVerts[vi] = new double[2];
	}

	/* find the first point of polygon j inside polygon i (assume clockwise orientation of vertices) */
	allInside = true;
	for( vj = 0; vj < poly->n; vj++ ) {
		if( IsInside( poly->verts[vj] ) && !IsInside( poly->verts[(poly->n+vj-1)%poly->n] ) ) {
			jOffset = vj;
			allInside = false;
			break;
		}
		if( vj == poly->n ) {
			jOffset = 0;
			break;
		}
	}

	/* polygon j is entirely inside polygon i */
	if( allInside ) {
		intPoly = new Polygon( poly->verts, poly->n, tris[0]->order );
		for( vi = 0; vi < MAX_POLY; vi++ ) {
			delete[] polyVerts[vi];
		}
		delete[] polyVerts;
		return intPoly;
	}

	/* find the first point of polygon i inside polygon j (assume clockwise orientation of vertices) */
	allInside = true;
	for( vi = 0; vi < n; vi++ ) {
		if( poly->IsInside( verts[vi] ) && !poly->IsInside( verts[(n+vi-1)%n] ) ) {
			iOffset = vi;
			allInside = false;
			break;
		}
		if( vi == n ) {
			iOffset = 0;
			break;
		}
	}

	/* polygon i is entirely inside polygon j */
	if( allInside ) {
		intPoly = new Polygon( verts, n, tris[0]->order );
		for( vi = 0; vi < MAX_POLY; vi++ ) {
			delete[] polyVerts[vi];
		}
		delete[] polyVerts;
		return intPoly;
	}

	/* no verts of either polygon inside the other */
	if( iOffset == -1 && jOffset == -1 ) {
		for( vi = 0; vi < MAX_POLY; vi++ ) {
			delete[] polyVerts[vi];
		}
		delete[] polyVerts;
		return NULL;
	}

	/* start adding polygon i's verts to the new polygon */
	if( iOffset > -1 ) {
		for( vi = 0; vi < n; vi++ ) {
			if( !poly->IsInside( verts[(vi+iOffset)%n] ) )
				break;

			pSize = AddPoint( polyVerts, pSize, verts[(vi+iOffset)%n] );
		}
	}

	/* add the ij intersection point to the new polygon */
	if( edges[(n+vi+iOffset-1)%n]->Intersection( poly->edges[(poly->n+jOffset-1)%poly->n], ijInt ) ) {
		pSize = AddPoint( polyVerts, pSize, ijInt );
	}

	/* start adding polygon j's verts to the new polygon */
	if( jOffset > -1 ) {
		for( vj = 0; vj < poly->n; vj++ ) {
			if( !IsInside( poly->verts[(vj+jOffset)%poly->n] ) )
				break;

			pSize = AddPoint( polyVerts, pSize, poly->verts[(vj+jOffset)%poly->n] );
		}
	}

	/* add the ji intersection point to the new polygon */
	if( poly->edges[(poly->n+vj+jOffset-1)%poly->n]->Intersection( edges[(n+iOffset-1)%n], jiInt ) ) {
		pSize = AddPoint( polyVerts, pSize, jiInt );
	}

	/* now instantiate the new polygon */
	if( pSize > 2 ) {
		intPoly = new Polygon( polyVerts, pSize, tris[0]->order );
	}

	for( vi = 0; vi < MAX_POLY; vi++ ) {
		delete[] polyVerts[vi];
	}
	delete[] polyVerts;

	return intPoly;
}

/* reqires that polygon vertices are consistently ordered clockwise */
bool Polygon::IsInside( double* pt ) {
	int         i;
	double      edgeCrossOrig, edgeCrossPt;
	Triangle*   tri;

	for( i = 0; i < n; i++ ) {
		/* is the point shared by both polygons? */
		if( (verts[i][0] - pt[0])*(verts[i][0] - pt[0]) + (verts[i][1] - pt[1])*(verts[i][1] - pt[1]) < R2_TOL ) {
			return true;
		}

		edgeCrossOrig = tris[i]->CrossProduct( 2 );
		tri = new Triangle( verts[i], verts[(i+1)%n], pt, 2 );
		edgeCrossPt = tri->CrossProduct( 2 );
		delete tri;

		if( edgeCrossOrig*edgeCrossPt < 0.0 ) {
			return false;
		}
	}

	return true;
}

void Polygon::GenOrigin() {
	int i;

	origin[0] = origin[1] = 0.0;

	for( i = 0; i < n; i++ ) {
		origin[0] += verts[i][0];
		origin[1] += verts[i][1];
	}
	origin[0] /= n;
	origin[1] /= n;
}

double Polygon::Area() {
	double  area = 0.0;
	int     tri_i;

	for( tri_i = 0; tri_i < n; tri_i++ ) {
		area += tris[tri_i]->area;
	}
	return area;
}

void Polygon::Print() {
	int vi;

	for( vi = 0; vi < n; vi++ ) {
		cout << "[" << verts[vi][0] << ", " << verts[vi][1] << "]" << endl;
	}
}
