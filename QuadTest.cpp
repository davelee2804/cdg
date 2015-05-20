#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Basis.h"
#include "Grid.h"
#include "Field.h"
#include "LinAlg.h"
#include "CFA.h"
#include "CDG.h"

using namespace std;

#define QUAD_ORDER 8
#define BASIS_ORDER 2

typedef double ( Func ) ( double* x );

double func( double* x ) {
	return cos( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
}

double Integrate( Field* field, Func* func ) {
    double      vol	    = 0.0;
    int         pi, ti, qi;
	Grid*		grid    = field->grid;
    Polygon*    poly;
    Triangle*   tri;

    for( pi = 0; pi < grid->nPolys; pi++ ) {
        poly = grid->polys[pi];
        for( ti = 0; ti < poly->n; ti++ ) {
            tri = poly->tris[ti];
            for( qi = 0; qi < tri->nq; qi++ ) {
                vol += tri->wi[qi]*func( tri->qi[qi] )*tri->area;
            }
        }
    }
    return vol;
}

int main() {
	int      nx     = 1;
	int      ny     = 1;
	int      i, j;
	Grid*    grid;
	Field*   field;
	CDG*     cdg;
	double   vol;
	double   ans    = 16.0/M_PI/M_PI;

	/* test convergence of errors with increased resolution for given quadrature order */
	cout << "testing quadrature, order: " << QUAD_ORDER << endl;
	for( i = 0; i < 8; i++ ) {
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );

		for( j = 0; j < grid->nPolys; j++ ) {
			field->basis[j]->ci[0] = func( grid->polys[j]->origin );
		}

		vol = Integrate( field, func );
		cout << fabs( 1.0 - vol/ans ) << endl;

		nx *= 2;
		ny *= 2;

		delete field;
		delete grid;
	}

	/* now test with the basis functions */
	nx = ny = 1;
	cout << "testing basis function setup and integration" << endl;
	for( i = 0; i < 8; i++ ) {
		grid  = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );
		cdg   = new CDG( field, NULL, NULL, NULL, NULL );

		cdg->InitBetaIJInv( func );
		vol = field->Integrate();
		cout << fabs( 1.0 - vol/ans ) << endl;

		nx *= 2;
		ny *= 2;

		delete cdg;
		delete field;
		delete grid;
	}

	return 1;
}
