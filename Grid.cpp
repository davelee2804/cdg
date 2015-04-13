#include <iostream>
#include <fstream>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Grid.h"

Grid::Grid( int _nx, int _ny, double _minx, double _miny, double _maxx, double _maxy, int _quadOrder, int _basisOrder, bool _internal ) {
	int 		i, j, k, l, ei = 0;
	double 		x, y;
	double** 	cellVerts;
	double**	points;
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

	nVerts = (nx+1)*(ny+1);
	nEdges = (nx+1)*ny + nx*(ny+1);
	nCells = nx*ny;

	verts = new double*[nVerts];
	edges = new Edge*[nEdges];
	cells = new Cell*[nCells];
	cellVerts = new double*[4];
	cellVerts[0] = new double[2];
	cellVerts[1] = new double[2];
	cellVerts[2] = new double[2];
	cellVerts[3] = new double[2];
	points = new double*[basisOrder*basisOrder];
	for( i = 0; i < basisOrder*basisOrder; i++ ) {
		points[i] = new double[2];
	}

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
			/* enter the verts for each cell clockwise */
			cellVerts[0][0] = verts[(j+0)*(nx+1)+(i+0)][0];
			cellVerts[0][1] = verts[(j+0)*(nx+1)+(i+0)][1];
			cellVerts[1][0] = verts[(j+1)*(nx+1)+(i+0)][0];
			cellVerts[1][1] = verts[(j+1)*(nx+1)+(i+0)][1];
			cellVerts[2][0] = verts[(j+1)*(nx+1)+(i+1)][0];
			cellVerts[2][1] = verts[(j+1)*(nx+1)+(i+1)][1];
			cellVerts[3][0] = verts[(j+0)*(nx+1)+(i+1)][0];
			cellVerts[3][1] = verts[(j+0)*(nx+1)+(i+1)][1];

			/* add the internal points for higher order representations */
			for( k = 0; k < basisOrder; k++ ) {
				for( l = 0; l < basisOrder; l++ ) {
					if( internal ) {
						points[k*basisOrder+l][0] = minx + i*dx + (0.5 + l)*dx/basisOrder;
						points[k*basisOrder+l][1] = miny + j*dy + (0.5 + k)*dy/basisOrder;
					}
					else {
						points[k*basisOrder+l][0] = minx + i*dx + l*dx/(basisOrder - 1);
						points[k*basisOrder+l][1] = miny + j*dy + k*dy/(basisOrder - 1);
					}
				}
			}

			/* enter the edges for each cell - order not important as searching routines use polygon verts and sub 
			   triangles, not edges. edge ordering is left/right/bottom/top. */
			cells[j*nx+i] = new Cell( cellVerts, 4, quadOrder, basisOrder*basisOrder, points );
		}
	}

	delete[] cellVerts[0];
	delete[] cellVerts[1];
	delete[] cellVerts[2];
	delete[] cellVerts[3];
	delete[] cellVerts;
	for( i = 0; i < basisOrder*basisOrder; i++ ) {
		delete[] points[i];
	}
	delete[] points;
}

Grid::~Grid() {
	int i;

	for( i = 0; i < nCells; i++ ) {
		delete cells[i];
	}
	for( i = 0; i < nEdges; i++ ) {
		delete edges[i];
	}
	for( i = 0; i < nVerts; i++ ) {
		delete[] verts[i];
	}
	delete[] cells;
	delete[] edges;
	delete[] verts;
}

Cell* Grid::FindCell( double* pt ) {
	int 	i 		= (pt[0] - minx)/dx;	
	int 	j 		= (pt[1] - miny)/dy;
	Cell* 	cell 	= cells[j*nx+i];

	/* sanity test */
	if( !cell->IsInside( pt ) ) {
		cerr << "ERROR in cell finding algorithm" << endl;
		abort();
	}
	return cell;
}

int Grid::GetCellIndex( double* pt ) {
	int 	i 		= (pt[0] - minx)/dx;	
	int 	j 		= (pt[1] - miny)/dy;

	return j*nx + i;
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

void Grid::GetEdgeCells( int ei, Cell* p1, Cell* p2 ) {
	int norm = ei/((nx+1)*ny);		// norm = 0: edge is normal to x, norm = 1: edge is normal to y
	int	dummy = ei - norm*(nx+1)*ny;
	int xi = ( norm == 0 ) ? dummy%(nx+1) : dummy%nx;
	int yj = ( norm == 0 ) ? dummy/(nx+1) : dummy/nx;

	/* polygons left and right */
	if( norm == 0 ) {
		/* disregard boundary cells for now */
		if( xi == 0 || xi == nx ) {	
			p1 = p2 = NULL;
			return;
		}
		p1 = cells[yj*nx+xi-1];
		p2 = cells[yj*nx+xi];
	}
	/* polygons bottom and top */
	else {
		/* disregard boundary cells for now */
		if( yj == 0 || yj == ny ) {	
			p1 = p2 = NULL;
			return;
		}
		p1 = cells[(yj-1)*nx+xi];
		p2 = cells[yj*nx+xi];
	}
}

bool Grid::GetEdgeCellInds( int ei, int* pinds ) {
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

void Grid::GetCellEdgeInds( int pi, int* einds ) {
	int xi 		= pi%nx;
	int yj 		= pi/nx;
	int shift 	= (nx+1)*ny;

	einds[0] = (yj+0)*(nx+1)+(xi+0);
	einds[1] = (yj+1)*(nx+0)+(xi+0) + shift;
	einds[2] = (yj+0)*(nx+1)+(xi+1);
	einds[3] = (yj+0)*(nx+0)+(xi+0) + shift;
}

void Grid::UpdateEdges() {
	int i, j, ei = 0;

	for( j = 0; j < ny; j++ ) {
		for( i = 0; i < nx+1; i++ ) {
			edges[ei]->v1[0] = verts[(j+0)*(nx+1)+(i+0)][0];
			edges[ei]->v1[1] = verts[(j+0)*(nx+1)+(i+0)][1];
			edges[ei]->v2[0] = verts[(j+1)*(nx+1)+(i+0)][0];
			edges[ei]->v2[1] = verts[(j+1)*(nx+1)+(i+0)][1];
			ei++;
		}
	}

	for( j = 0; j < ny+1; j++ ) {
		for( i = 0; i < nx; i++ ) {
			edges[ei]->v1[0] = verts[(j+0)*(nx+1)+(i+0)][0];
			edges[ei]->v1[1] = verts[(j+0)*(nx+1)+(i+0)][1];
			edges[ei]->v2[0] = verts[(j+0)*(nx+1)+(i+1)][0];
			edges[ei]->v2[1] = verts[(j+0)*(nx+1)+(i+1)][1];
			ei++;
		}
	}
}

void Grid::UpdateCells() {
	int i, j, pi = 0;

	for( j = 0; j < ny; j++ ) {
		for( i = 0; i < nx; i++ ) {
			cells[pi]->verts[0][0] = verts[(j+0)*(nx+1)+(i+0)][0];
			cells[pi]->verts[0][1] = verts[(j+0)*(nx+1)+(i+0)][1];
			cells[pi]->verts[1][0] = verts[(j+1)*(nx+1)+(i+0)][0];
			cells[pi]->verts[1][1] = verts[(j+1)*(nx+1)+(i+0)][1];
			cells[pi]->verts[2][0] = verts[(j+1)*(nx+1)+(i+1)][0];
			cells[pi]->verts[2][1] = verts[(j+1)*(nx+1)+(i+1)][1];
			cells[pi]->verts[3][0] = verts[(j+0)*(nx+1)+(i+1)][0];
			cells[pi]->verts[3][1] = verts[(j+0)*(nx+1)+(i+1)][1];
			pi++;
		}
	}
}

void Grid::Write( string fname ) {
	ofstream file;
	char filename[80];
	int i, j;

	sprintf( filename, "output/%s.x.txt", fname.c_str() );
	file.open( filename );
	for( i = 0; i < nx; i++ ) {
		for( j = 0; j < basisOrder; j++ ) {
			file << cells[i]->coords[j][0] << endl;
		}
	}
	file.close();

	sprintf( filename, "output/%s.y.txt", fname.c_str() );
	file.open( filename );
	for( i = 0; i < nCells; i += nx + 1 ) {
		for( j = 0; j < basisOrder*basisOrder; j += basisOrder ) {
			file << cells[i]->coords[j][1] << endl;
		}
	}
	file.close();
}
