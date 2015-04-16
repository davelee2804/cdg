#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"

using namespace std;

Field::Field( Grid* _grid ) {
	int i;

	grid = _grid;

	basis = new Basis*[grid->nCells];
	for( i = 0; i < grid->nCells; i++ ) {
		basis[i] = new Basis( grid->basisOrder, grid->cells[i]->origin, grid->cells[i] );
	}
}

Field::~Field() {
	int i;

	for( i = 0; i < grid->nCells; i++ ) {
		delete basis[i];
	}
	delete[] basis;
}

void Field::LinearInterp( double* x, double* v ) {
	int 	xi = (x[0] - grid->minx)/grid->dx;
	int 	yj = (x[1] - grid->miny)/grid->dy;
	double 	xl = x[0] - grid->minx - xi*grid->dx;
	double 	yl = x[1] - grid->miny - yj*grid->dy;
	int		i, j;
	double	xb[2], yb[2];

	if( grid->basisOrder != 2 ) {
		cerr << "Error! basis order is: " << grid->basisOrder << ", cannot perform linear interpolation" << endl;
		abort();
	}

	/* linear basis funcs */
	xb[0] = 1.0 - xl/grid->dx;
	xb[1] = xl/grid->dx;
	yb[0] = 1.0 - yl/grid->dy;
	yb[1] = yl/grid->dy;

	v[0] = 0.0;
	for( j = 0; j < 2; j++ ) {
		for( i = 0; i < 2; i++ ) {
			v[0] += basis[yj*grid->nx+xi]->ci[j*2+i]*xb[i]*yb[j];
		}
	}
}

/* integrate assuming constant values for each cell */
double Field::IntegrateConstant() {
	int i;
	double vol = 0;

	for( i = 0; i < grid->nCells; i++ ) {
		vol += basis[i]->ci[0]*grid->dx*grid->dy;
	}
	return vol;
}

double Field::Integrate() {
	int i, j, k;
	double val, vol = 0;
	Cell* cell;
	Triangle* tri;

	for( i = 0; i < grid->nCells; i++ ) {
		cell = grid->cells[i];
		for( j = 0; j < cell->n; j++ ) {
			tri = cell->tris[j];
			for( k = 0; k < tri->nQuadPts; k++ ) {
				val = basis[i]->EvalFull( tri->qi[k] );
				vol += val*tri->wi[k]*tri->Area();
			}
		}
	}
	return vol;
}

void Field::Copy( Field* field ) {
	int i, j;

	for( i = 0; i < grid->nCells; i++ ) {
		for( j = 0; j < grid->cells[i]->nc; j++ ) {
			basis[i]->ci[j] = field->basis[i]->ci[j];
		}
	}
}

void Field::Write( string fname, int tstep ) {
	ofstream file;
	char filename[80];
	int i, j, k, l;
	Cell* cell;
	Basis* basis_i;

	sprintf( filename, "output/%s.%.4u.txt", fname.c_str(), tstep );
	file.open( filename );
	for( i = 0; i < grid->ny; i++ ) {
		for( j = 0; j < grid->basisOrder; j++ ) {
			for( k = 0; k < grid->nx; k++ ) {
				for( l = 0; l < grid->basisOrder; l++ ) {
					cell = grid->cells[i*grid->nx+k];
					basis_i = basis[i*grid->nx+k];
					file << basis_i->EvalFull( cell->coords[j*grid->basisOrder+l] ) << endl;
				}
			}
		}
		file << endl;
	}
	file.close();
}
