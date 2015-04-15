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
#include "CFA.h"
#include "CDG.h"

using namespace std;

#define QUAD_ORDER 2
#define BASIS_ORDER 2

double func( double* x ) {
	return cos( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
}

/* test the basis function representation. since the grid is cartesian, use a first order basis:
 * 		c_0 + c_1.x + c_2.y + c_3.xy, 
 * located at the first order quadrature points for the quadralateral */
int main() {
	int			nx			= 1;
	int			ny			= 1;
	int			i, j, k;
	Grid*		grid;
	Field*		field;
	Cell* 		cell;
	Basis*		basis;
	double		ans			= 16.0/M_PI/M_PI;
	double		vol;
	CDG*		cdg;
	Field*		velx;
	Field*		vely;
	double*		beta_ij_2;
	double*		cj;

	cout << "testing the basis function matrix inverse..." << endl;
	for( i = 0; i < 8; i++ ) {
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );
		velx = new Field( grid );
		vely = new Field( grid );

		for( j = 0; j < grid->nCells; j++ ) {
			for( k = 0; k < grid->cells[j]->nc; k++ ) {
				field->basis[j]->ci[k] = func( grid->cells[j]->coords[k] );
			}
		}
		cdg = new CDG( field, velx, vely );

		beta_ij_2 = new double[grid->cells[0]->nc*grid->cells[0]->nc];
		cj = new double[grid->cells[0]->nc];

		for( j = 0; j < grid->nCells; j++ ) {
			cell = grid->cells[j];
			basis = field->basis[j];

			MatInv( cdg->betaInv_ij[j], beta_ij_2, cell->nc );
			AXEB( beta_ij_2, basis->ci, cj, cell->nc );
			for( k = 0; k < cell->nc; k++ ) {
				basis->ci[k] = cj[k];
			}
		}

		vol = field->Integrate();

		cout << fabs( 1.0 - vol/ans ) << endl;

		delete[] beta_ij_2;
		delete[] cj;
		delete cdg;
		delete grid;
		delete field;
		delete velx;
		delete vely;

		nx *= 2;
		ny *= 2;
	}

	return 1;
}
