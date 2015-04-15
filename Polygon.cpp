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

	origin[0] = 0.0;
	origin[1] = 0.0;
	for( i = 0; i < n; i++ ) {
		verts[i] = new double[2];
		verts[i][0] = _verts[i][0];
		verts[i][1] = _verts[i][1];
		origin[0] += verts[i][0];
		origin[1] += verts[i][1];
	}
	origin[0] /= n;
	origin[1] /= n;

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

#define R2_TOL 1.0e-6

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
	int 		vi, vj, iOffset = -1, jOffset = -1, pSize = 0;
	double		**polyVerts;
	double		ijInt[2], jiInt[2];
	Polygon*	intPoly	= NULL;
	bool		allInside;

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
	int 		i;
	double		edgeCrossOrig, edgeCrossPt;
	Triangle*	tri;

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

double Polygon::Area() {
	double 	area 	= 0.0;
	int 	tri_i;

	for( tri_i = 0; tri_i < n; tri_i++ ) {
		area += tris[tri_i]->Area();
	}
	return area;
}

void Polygon::Print() {
	int vi;

	for( vi = 0; vi < n; vi++ ) {
		cout << "[" << verts[vi][0] << ", " << verts[vi][1] << "]" << endl;
	}
}

#define I(p,q,N) \
        p*N + q
 
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void MatInv( double* A, double* Ainv, int n ) {
	int *indxc, *indxr, *ipiv;
	int i, j, k, l, irow = 0, icol = 0, ll;
	double big, dum, pivinv, temp;
 
	indxc = new int[n]; indxr = new int[n]; ipiv  = new int[n];
 
	for( i = 0; i < n*n; i++ ) { Ainv[i] = A[i]; }
	for( j = 0; j< n; j++ ) { ipiv[j] = 0; }
	for( i = 0; i < n; i++ ) {
		big = 0.0;
		for( j = 0; j < n; j++ ) {
			if( ipiv[j] != 1 ) {
				for( k = 0; k < n; k++ ) {
					if( ipiv[k] == 0 ) {
						if( fabs(Ainv[I(j,k,n)]) >= big ) {
							big = fabs(Ainv[I(j,k,n)]);
							irow = j;
							icol = k;
						}
					}
					else if( ipiv[k] > 1 ) { cerr << "Matrix inverse error! - singular matrix (1)\n"; }
				}
			}
		}
		++(ipiv[icol]);
		if( irow != icol ) {
			for( l = 0; l < n; l++ ) {
        			SWAP( Ainv[I(irow,l,n)], Ainv[I(icol,l,n)] );
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if( fabs(Ainv[I(icol,icol,n)]) < 1.0e-12 ) { cerr << "Matrix inverse error! - singular matrix (2)\n"; }
		pivinv = 1.0/Ainv[I(icol,icol,n)];
		Ainv[I(icol,icol,n)] = 1.0;
		for( l = 0; l < n; l++ ) { Ainv[I(icol,l,n)] *= pivinv; }
		for( ll = 0; ll < n; ll++ ) {
			if( ll != icol ) {
				dum = Ainv[I(ll,icol,n)];
				Ainv[I(ll,icol,n)] = 0.0;
				for( l = 0; l < n; l++ ) { Ainv[I(ll,l,n)] -= Ainv[I(icol,l,n)]*dum; }
			}
		}
	}
	for( l = n-1; l >= 0; l-- ) {
		if( indxr[l] != indxc[l] ) {
			for( k = 0; k < n; k++ ) {
				SWAP( Ainv[I(k,indxr[l],n)], Ainv[I(k,indxc[l],n)] );
			}
		}
	}
	delete[] indxc; delete[] indxr; delete[] ipiv;
}

/* b_i = A_{ij}x_j */
void AXEB( double* A, double* x, double* b, int n ) {
	int i, j;
 
	for( i = 0; i < n; i++ ) {
		b[i] = 0.0;
		for( j = 0; j < n; j++ ) {
			b[i] += A[I(i,j,n)]*x[j];
		}
	}
}
