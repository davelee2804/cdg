#include <iostream>
#include <fstream>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Grid.h"

Grid::Grid( int _nx, int _ny, double _minx, double _miny, double _maxx, double _maxy, int _quadOrder, int _basisOrder, bool _internal ) {
	int 		i, j, ei = 0;
	double 		x, y;
	double** 	polyVerts;
	double		*p1, *p2;

	nx = _nx;
	ny = _ny;
	minx = _minx;
	miny = _miny;
	maxx = _maxx;
	maxy = _maxy;
	quadOrder = _quadOrder;
	basisOrder = _basisOrder;
	internal = _internal;

	dx = (maxx - minx)/nx;
	dy = (maxy - miny)/ny;
	dxInv = 1.0/dx;
	dyInv = 1.0/dy;

	nVerts = (nx+1)*(ny+1);
	nEdges = (nx+1)*ny + nx*(ny+1);
	nPolys = nx*ny;

	verts = new double*[nVerts];
	edges = new Edge*[nEdges];
	polys = new Polygon*[nPolys];
	polyVerts = new double*[4];
	polyVerts[0] = new double[2];
	polyVerts[1] = new double[2];
	polyVerts[2] = new double[2];
	polyVerts[3] = new double[2];

	/* generate the verts */
	for( j = 0; j < ny+1; j++ ) {
		y = miny + j*dy;
		for( i = 0; i < nx+1; i++ ) {
			x = minx + i*dx;
			verts[j*(nx+1)+i] = new double[2];
			verts[j*(nx+1)+i][0] = x;
			verts[j*(nx+1)+i][1] = y;
		}
	}

	/* generate the edges - do those normal to x first */
	for( j = 0; j < ny; j++ ) {
		for( i = 0; i < nx+1; i++ ) {
			p1 = verts[(j+0)*(nx+1)+(i+0)];
			p2 = verts[(j+1)*(nx+1)+(i+0)];
			edges[ei++] = new Edge( p1, p2 );
		}
	}
	/* then do the ones normal to y */
	for( j = 0; j < ny+1; j++ ) {
		for( i = 0; i < nx; i++ ) {
			p1 = verts[(j+0)*(nx+1)+(i+0)];
			p2 = verts[(j+0)*(nx+1)+(i+1)];
			edges[ei++] = new Edge( p1, p2 );
		}
	}

	/* generate the polygons */
	for( j = 0; j < ny; j++ ) {
		for( i = 0; i < nx; i++ ) {
			/* enter the verts for each poly clockwise */
			polyVerts[0][0] = verts[(j+0)*(nx+1)+(i+0)][0];
			polyVerts[0][1] = verts[(j+0)*(nx+1)+(i+0)][1];
			polyVerts[1][0] = verts[(j+1)*(nx+1)+(i+0)][0];
			polyVerts[1][1] = verts[(j+1)*(nx+1)+(i+0)][1];
			polyVerts[2][0] = verts[(j+1)*(nx+1)+(i+1)][0];
			polyVerts[2][1] = verts[(j+1)*(nx+1)+(i+1)][1];
			polyVerts[3][0] = verts[(j+0)*(nx+1)+(i+1)][0];
			polyVerts[3][1] = verts[(j+0)*(nx+1)+(i+1)][1];

			/* enter the edges for each poly - order not important as searching routines use polygon verts and sub 
			   triangles, not edges. edge ordering is left/right/bottom/top. */
			polys[j*nx+i] = new Polygon( polyVerts, 4, quadOrder );
		}
	}

	delete[] polyVerts[0];
	delete[] polyVerts[1];
	delete[] polyVerts[2];
	delete[] polyVerts[3];
	delete[] polyVerts;
}

Grid::~Grid() {
	int i;

	for( i = 0; i < nPolys; i++ ) {
		delete polys[i];
	}
	for( i = 0; i < nEdges; i++ ) {
		delete edges[i];
	}
	for( i = 0; i < nVerts; i++ ) {
		delete[] verts[i];
	}
	delete[] polys;
	delete[] edges;
	delete[] verts;
}

int Grid::GetPolyIndex( double* pt ) {
	int 	i 		= (pt[0] - minx)/dx;	
	int 	j 		= (pt[1] - miny)/dy;
	int		ind		= j*nx + i;

	if( polys[ind]->IsInside( pt ) ) {
		return ind;
	}

	if( i > 0      && j > 0      && polys[(j-1)*nx+(i-1)]->IsInside( pt ) ) {
		return (j-1)*nx + (i-1);
	}
	if( i > 0      &&               polys[(j+0)*nx+(i-1)]->IsInside( pt ) ) {
		return (j+0)*nx + (i-1);
	}
	if( i > 0      && j < ny - 1 && polys[(j+1)*nx+(i-1)]->IsInside( pt ) ) {
		return (j+1)*nx + (i-1);
	}
	if(               j < ny - 1 && polys[(j+1)*nx+(i+0)]->IsInside( pt ) ) {
		return (j+1)*nx + (i+0);
	}
	if( i < nx - 1 && j < ny - 1 && polys[(j+1)*nx+(i+1)]->IsInside( pt ) ) {
		return (j+1)*nx + (i+1);
	}
	if( i < nx - 1 &&               polys[(j+0)*nx+(i+1)]->IsInside( pt ) ) {
		return (j+0)*nx + (i+1);
	}
	if( i < nx - 1 && j > 0      && polys[(j-1)*nx+(i+1)]->IsInside( pt ) ) {
		return (j-1)*nx + (i+1);
	}
	if(               j > 0      && polys[(j-1)*nx+(i+0)]->IsInside( pt ) ) {
		return (j-1)*nx + (i+0);
	}

	return ind;
}

/* norm = 0: edge is normal to x, norm = 1: edge is normal to y */
int Grid::EdgeCoordToIndex( int norm, int xi, int yj ) {
	int offset = norm*(nx+1)*ny;
	int	index;

	if( norm == 0 ) {
		index = yj*(nx+1) + xi;
	}
	else if( norm == 1 ) {
		index = yj*nx + xi;
	}
	else {
		cerr << "invalid normal index in Grid::EdgeCoordToIndex, must 0 or 1" << endl;
		abort();
	}
	return offset + index;
}

void Grid::EdgeIndexToCoord( int ei, int* norm, int* xi, int* yj ) {
	int	dummy;

	*norm = ei/((nx+1)*ny);		// norm = 0: edge is normal to x, norm = 1: edge is normal to y
	
	dummy = ei - (*norm)*(nx+1)*ny;
	if( *norm == 0 ) {
		*xi = dummy%(nx+1);
		*yj = dummy/(nx+1);
	}
	else {
		*xi = dummy%nx;
		*yj = dummy/nx;
	}
}

void Grid::GetEdgePolys( int ei, Polygon* p1, Polygon* p2 ) {
	int norm = ei/((nx+1)*ny);		// norm = 0: edge is normal to x, norm = 1: edge is normal to y
	int	dummy = ei - norm*(nx+1)*ny;
	int xi = ( norm == 0 ) ? dummy%(nx+1) : dummy%nx;
	int yj = ( norm == 0 ) ? dummy/(nx+1) : dummy/nx;

	/* polygons left and right */
	if( norm == 0 ) {
		/* disregard boundary polys for now */
		if( xi == 0 || xi == nx ) {	
			p1 = p2 = NULL;
			return;
		}
		p1 = polys[yj*nx+xi-1];
		p2 = polys[yj*nx+xi];
	}
	/* polygons bottom and top */
	else {
		/* disregard boundary polys for now */
		if( yj == 0 || yj == ny ) {	
			p1 = p2 = NULL;
			return;
		}
		p1 = polys[(yj-1)*nx+xi];
		p2 = polys[yj*nx+xi];
	}
}

bool Grid::GetEdgePolyInds( int ei, int* pinds ) {
	int norm, xi, yj;

	EdgeIndexToCoord( ei, &norm, &xi, &yj );

	/* ignore boundaries for now */
	if( ( xi == 0 ) || ( yj == 0 ) || ( norm == 0 && xi == nx ) || ( norm == 1 && yj == ny ) ) {
		return false;
	}
	/* ignore edges incident on boundary also */
	if( ( norm == 0 && yj == 0 ) || ( norm == 0 && yj == ny - 1 ) || ( norm == 1 && xi == 0 ) || ( norm == 1 && yj == nx - 1 ) ) {
		return false;
	}

	if( norm == 0 ) {
		pinds[0] = (yj-1)*nx+(xi-1);
		pinds[1] = (yj-1)*nx+(xi+0);
		pinds[2] = (yj+0)*nx+(xi-1);
		pinds[3] = (yj+0)*nx+(xi+0);
		pinds[4] = (yj+1)*nx+(xi-1);
		pinds[5] = (yj+1)*nx+(xi+0);
	}
	else {
		pinds[0] = (yj-1)*nx+(xi-1);
		pinds[1] = (yj-1)*nx+(xi+0);
		pinds[2] = (yj-1)*nx+(xi+1);
		pinds[3] = (yj+0)*nx+(xi-1);
		pinds[4] = (yj+0)*nx+(xi+0);
		pinds[5] = (yj+0)*nx+(xi+1);
	}

	return true;
}

void Grid::GetPolyEdgeInds( int pi, int* einds ) {
	int xi 		= pi%nx;
	int yj 		= pi/nx;
	int shift 	= (nx+1)*ny;

	einds[0] = (yj+0)*(nx+1)+(xi+0);
	einds[1] = (yj+1)*(nx+0)+(xi+0) + shift;
	einds[2] = (yj+0)*(nx+1)+(xi+1);
	einds[3] = (yj+0)*(nx+0)+(xi+0) + shift;
}

/* no chech for bounndary vertices yet */
void Grid::GetVertPolyInds( int vi, int* cinds ) {
	int xi = vi%(nx+1);
	int yj = vi/(nx+1);

	cinds[0] = (yj-1)*nx + (xi-1);
	cinds[1] = (yj-1)*nx + (xi+0);
	cinds[2] = (yj+0)*nx + (xi-1);
	cinds[3] = (yj+0)*nx + (xi+0);
}

/* return in clockwise orientation */
void Grid::GetPolyVertInds( int ci, int* vinds ) {
	int xi = ci%nx;
	int yj = ci/nx;

	vinds[0] = (yj+0)*(nx+1) + (xi+0);
	vinds[1] = (yj+1)*(nx+1) + (xi+0);
	vinds[2] = (yj+1)*(nx+1) + (xi+1);
	vinds[3] = (yj+0)*(nx+1) + (xi+1);
}

void Grid::UpdateEdges() {
	int i, j, ei = 0;

	for( j = 0; j < ny; j++ ) {
		for( i = 0; i < nx+1; i++ ) {
			edges[ei]->v1[0] = verts[(j+0)*(nx+1)+(i+0)][0];
			edges[ei]->v1[1] = verts[(j+0)*(nx+1)+(i+0)][1];
			edges[ei]->v2[0] = verts[(j+1)*(nx+1)+(i+0)][0];
			edges[ei]->v2[1] = verts[(j+1)*(nx+1)+(i+0)][1];
			edges[ei]->Init();
			ei++;
		}
	}

	for( j = 0; j < ny+1; j++ ) {
		for( i = 0; i < nx; i++ ) {
			edges[ei]->v1[0] = verts[(j+0)*(nx+1)+(i+0)][0];
			edges[ei]->v1[1] = verts[(j+0)*(nx+1)+(i+0)][1];
			edges[ei]->v2[0] = verts[(j+0)*(nx+1)+(i+1)][0];
			edges[ei]->v2[1] = verts[(j+0)*(nx+1)+(i+1)][1];
			edges[ei]->Init();
			ei++;
		}
	}
}

void Grid::UpdatePolys() {
	Polygon* poly;
	int i, j, k, pi = 0;

	for( j = 0; j < ny; j++ ) {
		for( i = 0; i < nx; i++ ) {
			poly = polys[pi++];
			poly->verts[0][0] = verts[(j+0)*(nx+1)+(i+0)][0];
			poly->verts[0][1] = verts[(j+0)*(nx+1)+(i+0)][1];
			poly->verts[1][0] = verts[(j+1)*(nx+1)+(i+0)][0];
			poly->verts[1][1] = verts[(j+1)*(nx+1)+(i+0)][1];
			poly->verts[2][0] = verts[(j+1)*(nx+1)+(i+1)][0];
			poly->verts[2][1] = verts[(j+1)*(nx+1)+(i+1)][1];
			poly->verts[3][0] = verts[(j+0)*(nx+1)+(i+1)][0];
			poly->verts[3][1] = verts[(j+0)*(nx+1)+(i+1)][1];
			poly->GenOrigin();
			for( k = 0; k < poly->n; k++ ) {
				/* update the polygon edges */
				poly->edges[k]->v1[0] = poly->verts[k][0];
				poly->edges[k]->v1[1] = poly->verts[k][1];
				poly->edges[k]->v2[0] = poly->verts[(k+1)%poly->n][0];
				poly->edges[k]->v2[1] = poly->verts[(k+1)%poly->n][1];
				poly->edges[k]->Init();
				/* update the polygon triangles */
				poly->tris[k]->a[0] = poly->verts[k][0];
				poly->tris[k]->a[1] = poly->verts[k][1];
				poly->tris[k]->b[0] = poly->verts[(k+1)%poly->n][0];
				poly->tris[k]->b[1] = poly->verts[(k+1)%poly->n][1];
				poly->tris[k]->c[0] = poly->origin[0];
				poly->tris[k]->c[1] = poly->origin[1];
				poly->tris[k]->Init();
			}
		}
	}
}

void Grid::UpdateTris() {
	int 		i, j;
	Polygon* 	poly;
	Triangle*	tri;

	for( i = 0; i < nPolys; i++ ) {
		poly = polys[i];
		for( j = 0; j < poly->n; j++ ) {
			tri = poly->tris[j];
			tri->a[0] = poly->verts[j][0];
			tri->a[1] = poly->verts[j][1];
			tri->b[0] = poly->verts[(j+1)%poly->n][0];
			tri->b[1] = poly->verts[(j+1)%poly->n][1];
			tri->c[0] = poly->origin[0];
			tri->c[1] = poly->origin[1];
			tri->Init();
			/* update the triangle edges */
			tri->ab->v1[0] = tri->a[0];
			tri->ab->v1[1] = tri->a[1];
			tri->ab->v2[0] = tri->b[0];
			tri->ab->v2[1] = tri->b[1];

			tri->bc->v1[0] = tri->b[0];
			tri->bc->v1[1] = tri->b[1];
			tri->bc->v2[0] = tri->c[0];
			tri->bc->v2[1] = tri->c[1];

			tri->ca->v1[0] = tri->c[0];
			tri->ca->v1[1] = tri->c[1];
			tri->ca->v2[0] = tri->a[0];
			tri->ca->v2[1] = tri->a[1];

			tri->Init();
		}
	}
}

void Grid::Write( string fname, int n ) {
	ofstream 	file;
	char 		filename[80];
	int 		i, j;
	double*		xCoords	= new double[n];
	double*		yCoords	= new double[n];

	if( internal ) {
		for( i = 0; i < n; i++ ) {
			xCoords[i] = (0.5 + i)*dx/n;
			yCoords[i] = (0.5 + i)*dy/n;
		}
	}
	else {
		for( i = 0; i < n; i++ ) {
			xCoords[i] = i*dx/(n-1);
			yCoords[i] = i*dy/(n-1);
		}
	}

	sprintf( filename, "output/%s.x.txt", fname.c_str() );
	file.open( filename );
	for( i = 0; i < nx; i++ ) {
		for( j = 0; j < n; j++ ) {
			file << minx + i*dx + xCoords[j] << endl;
		}
	}
	file.close();

	sprintf( filename, "output/%s.y.txt", fname.c_str() );
	file.open( filename );
	for( i = 0; i < ny; i++ ) {
		for( j = 0; j < n; j ++ ) {
			file << miny + i*dy + yCoords[j] << endl;
		}
	}
	file.close();

	delete[] xCoords;
	delete[] yCoords;
}

void Grid::WriteTris( string fname ) {
	ofstream   file;
	char       filename[80];
	Polygon*   poly;
	Triangle*  tri;
	int        i, j, k;

	sprintf( filename, "output/%s.txt", fname.c_str() );
	file.open( filename );
	for( i = 0; i < nPolys; i++ ) {
		poly = polys[i];
		for( j = 0; j < poly->n; j++ ) {
			tri = poly->tris[j];
			for( k = 0; k < tri->nq; k++ ) {
				file << i << "\t" << tri->qi[k][0] << "\t" << tri->qi[k][1] << endl;
			}
		}
	}
	file.close();
}
