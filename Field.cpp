#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"

using namespace std;

Field::Field( Grid* _grid ) {
	int i;

	grid = _grid;

	basis = new Basis*[grid->nPolys];
	for( i = 0; i < grid->nPolys; i++ ) {
		basis[i] = new Basis( grid->polys[i], grid->basisOrder, grid->polys[i]->origin, grid->dx, grid->dy );
	}
}

Field::~Field() {
	int i;

	for( i = 0; i < grid->nPolys; i++ ) {
		delete basis[i];
	}
	delete[] basis;
}

double Field::EvalAtCoord( double* x ) {
	int pi = grid->GetPolyIndex( x );

	return basis[pi]->EvalFull( x ); 
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

/* integrate assuming constant values for each poly */
double Field::IntegrateConstant() {
	int i;
	double vol = 0.0;

	for( i = 0; i < grid->nPolys; i++ ) {
		vol += basis[i]->ci[0]*grid->dx*grid->dy;
	}
	return vol;
}

double Field::Integrate() {
	int i, j, k;
	double val, vol = 0.0;
	Polygon* poly;
	Triangle* tri;

	for( i = 0; i < grid->nPolys; i++ ) {
		poly = grid->polys[i];
		for( j = 0; j < poly->n; j++ ) {
			tri = poly->tris[j];
			for( k = 0; k < tri->nQuadPts; k++ ) {
				val = basis[i]->EvalFull( tri->qi[k] );
				vol += val*tri->wi[k]*tri->Area();
			}
		}
	}
	return vol;
}

double Field::L1Error( Func* analytic, bool doNorm ) {
	Polygon*	poly;
	Triangle*	tri;
	double		error = 0.0, norm = 0.0;
	double		weight, val, ana;
	int 		i, j, k;

	for( i = 0; i < grid->nPolys; i++ ) {
		poly = grid->polys[i];
		for( j = 0; j < poly->n; j++ ) {
			tri = poly->tris[j];
			for( k = 0; k < tri->nQuadPts; k++ ) {
				val = basis[i]->EvalFull( tri->qi[k] );
				ana = analytic( tri->qi[k] );
				weight = tri->wi[k]*tri->Area();
				error += weight*fabs( val - ana );
				norm += weight*fabs( ana );
			}
		}
	}

	if( !doNorm ) {
		norm = 1.0;
	}

	return error / norm;
}

double Field::L2Error( Func* analytic, bool doNorm ) {
	Polygon*	poly;
	Triangle*	tri;
	double		errorSq = 0.0, normSq = 0.0;
	double		weight, val, ana;
	int 		i, j, k;

	for( i = 0; i < grid->nPolys; i++ ) {
		poly = grid->polys[i];
		for( j = 0; j < poly->n; j++ ) {
			tri = poly->tris[j];
			for( k = 0; k < tri->nQuadPts; k++ ) {
				val = basis[i]->EvalFull( tri->qi[k] );
				ana = analytic( tri->qi[k] );
				weight = tri->wi[k]*tri->Area();
				errorSq += weight*( val - ana )*( val - ana );
				normSq += weight*ana*ana;
			}
		}
	}

	if( !doNorm ) {
		normSq = 1.0;
	}

	return sqrt( errorSq / normSq );
}

void Field::Copy( Field* field ) {
	int i, j;

	for( i = 0; i < grid->nPolys; i++ ) {
		for( j = 0; j < basis[i]->nFuncs; j++ ) {
			basis[i]->ci[j] = field->basis[i]->ci[j];
		}
	}
}

void Field::UpdateBasis() {
	int i;

	for( i = 0; i < grid->nPolys; i++ ) {
		basis[i]->aInv = 1.0/grid->polys[i]->Area();
	}
}

void Field::Write( string fname, int tstep, int n ) {
	ofstream 	file;
	char 		filename[80];
	int 		i, j, k, l;
	Basis* 		basis_i;
	double		point[2];

	sprintf( filename, "output/%s.%.4u.txt", fname.c_str(), tstep );
	file.open( filename );
	for( i = 0; i < grid->ny; i++ ) {
		for( j = 0; j < n; j++ ) {
			for( k = 0; k < grid->nx; k++ ) {
				for( l = 0; l < n; l++ ) {
					basis_i = basis[i*grid->nx+k];
					if( grid->internal ) {
						point[0] = grid->minx + k*grid->dx + (0.5 + l)*grid->dx/n;
						point[1] = grid->miny + i*grid->dy + (0.5 + j)*grid->dy/n;
					}
					else {
						point[0] = grid->minx + k*grid->dx + l*grid->dx/(n-1);
						point[1] = grid->miny + i*grid->dy + j*grid->dy/(n-1);
					}
					file << basis_i->EvalFull( point ) << endl;
				}
			}
		}
		file << endl;
	}
	file.close();
}

void Field::WriteBasis( string fname, int tstep ) {
	ofstream 	file;
	char 		filename[80];
	int 		i, j;

	sprintf( filename, "output/%s_basis.%.4u.txt", fname.c_str(), tstep );
	file.open( filename );
	for( i = 0; i < grid->nPolys; i++ ) {
		for( j = 0; j < basis[i]->nFuncs; j++ ) {
			file << basis[i]->ci[j] << "\t";
		}
		file << endl;
	}
	file.close();
}

void Field::ReadBasis( string fname ) {
	ifstream 	file;
	int 		i, j;

	file.open( fname.c_str() );
	for( i = 0; i < grid->nPolys; i++ ) {
		for( j = 0; j < basis[i]->nFuncs; j++ ) {
			file >> basis[i]->ci[j];
		}
	}
	file.close();
}
