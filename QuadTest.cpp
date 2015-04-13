#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Basis.h"
#include "Grid.h"
#include "Field.h"

using namespace std;

#define QUAD_ORDER 2
#define BASIS_ORDER 1

typedef double ( Func ) ( double* x );

double func( double* x ) {
	return cos( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
}

double Integrate( Field* field, Func* func ) {
    double      volg    = 0.0;
    double      volp;
    double      volt;
    int         pi, ti, qi;
	Grid*		grid	= field->grid;
    Cell*    	cell;
    Triangle*   tri;

    for( pi = 0; pi < grid->nCells; pi++ ) {
        volp = 0.0;
        cell = grid->cells[pi];
        for( ti = 0; ti < cell->n; ti++ ) {
            volt = 0.0;
            tri = cell->tris[ti];
            for( qi = 0; qi < tri->nQuadPts; qi++ ) {
                volt += tri->wi[qi]*func( tri->qi[qi] )*tri->Area();
            }
            volp += volt;
        }

        volg += volp;
    }
    return volg;
}

int main() {
	int			nx			= 1;
	int			ny			= 1;
	int			i, j;
	Grid*		grid;
	Field*		field;
	double		vol;
	double		ans			= 16.0/M_PI/M_PI;

	for( i = 0; i < 8; i++ ) {
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );

		for( j = 0; j < grid->nCells; j++ ) {
			field->basis[j]->ci[0] = func( grid->cells[j]->origin );
		}

		vol = Integrate( field, func );
		cout << fabs( 1.0 - vol/ans ) << endl;

		nx *= 2;
		ny *= 2;

		delete grid;
		delete field;
	}

	return 1;
}
