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

#define QUAD_ORDER 1
#define BASIS_ORDER 2

double func( double* x ) {
	return cos( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
}

void GenBetaIPointJ( Basis* basis, double** pts, double* beta_ij ) {
    int     beta_i, beta_j;
    int     order   = basis->order;
    int     nBasis  = order*order;

    for( beta_j = 0; beta_j < nBasis; beta_j++ ) {
        for( beta_i = 0; beta_i < nBasis; beta_i++ ) {
            beta_ij[beta_j*nBasis+beta_i] = 0.0;
        }
    }

    for( beta_j = 0; beta_j < nBasis; beta_j++ ) {
        for( beta_i = 0; beta_i < nBasis; beta_i++ ) {
            beta_ij[beta_j*nBasis+beta_i] = basis->EvalIJ( pts[beta_j], beta_i, pts[beta_i] );
        }
    }
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
	double		**pts;
	int			nBasis		= BASIS_ORDER*BASIS_ORDER;
	double		fj[nBasis];
	double		beta_ij[nBasis*nBasis];
	double		betaInv_ij[nBasis*nBasis];
	CDG*		cdg;
	Field*		velx;
	Field*		vely;
	double*		beta_ij_2;
	//double		ci[4];
	double		*cj1, *cj2;
	double		weight;

	pts = new double*[4];
	for( k = 0; k < 4; k++ ) {
		pts[k] = new double[2];
	}

	cout << "testing the basis functions..." << endl;

	for( i = 0; i < 10; i++ ) {
		vol = 0.0;
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );

		for( j = 0; j < grid->nCells; j++ ) {
			cell = grid->cells[j];
			basis = field->basis[j];

			for( k = 0; k < 4; k++ ) {
				basis->ci[k] = func( cell->tris[k]->qi[0] );
				//ci[k] = func( cell->tris[k]->qi[0] );
				pts[k][0] = cell->tris[k]->qi[0][0];
				pts[k][1] = cell->tris[k]->qi[0][1];
				fj[k] = func( cell->tris[k]->qi[0] );
			}

			GenBetaIPointJ( basis, pts, beta_ij );
			MatInv( beta_ij, betaInv_ij, nBasis );
			AXEB( betaInv_ij, fj, basis->ci, nBasis );
			//AXEB( betaInv_ij, fj, ci, nBasis );

			for( k = 0; k < 4; k++ ) {
				vol += basis->EvalFull( cell->tris[k]->qi[0], pts )*cell->tris[k]->wi[0]*cell->tris[k]->Area();
				//vol += ci[k]*cell->tris[k]->wi[0]*cell->tris[k]->Area();
			}
		}

		cout << fabs( 1.0 - vol/ans ) << endl;

		nx *= 2;
		ny *= 2;

		delete grid;
		delete field;
	}

	/* now test the cdg beta_ij matrix construction */
	cout << "testing the basis function matrix inverse..." << endl;
	nx = 1;
	ny = 1;
	for( i = 0; i < 10; i++ ) {
		vol = 0.0;
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, 2, 2, true );
		field = new Field( grid );
		velx = new Field( grid );
		vely = new Field( grid );

		for( j = 0; j < grid->nCells; j++ ) {
			for( k = 0; k < grid->cells[0]->nc; k++ ) {
				field->basis[j]->ci[k] = func( grid->cells[j]->coords[k] );
			}
		}

		cdg = new CDG( field, velx, vely );

		beta_ij_2 = new double[grid->cells[0]->nc*grid->cells[0]->nc];
		cj1 = new double[grid->cells[0]->nc];
		cj2 = new double[grid->cells[0]->nc];

		for( j = 0; j < grid->nCells; j++ ) {
			cell = grid->cells[j];
			basis = field->basis[j];

			for( k = 0; k < cell->nc; k++ ) {
				cj1[k] = func( cell->coords[k] );
			}

			MatInv( cdg->betaInv_ij[j], beta_ij_2, cell->nc );
			AXEB( beta_ij_2, basis->ci, cj2, cell->nc );

			for( k = 0; k < cell->nc; k++ ) {
				/* basis functions are evaluated at the quadrature points, however fields are represented at the cell points... */
				weight = grid->dx*grid->dy/basis->order/basis->order;
				vol += weight*cj2[k];
			}
		}

		cout << fabs( 1.0 - vol/ans ) << endl;

		delete[] beta_ij_2;
		delete[] cj1;
		delete[] cj2;
		delete cdg;
		delete grid;
		delete field;
		delete velx;
		delete vely;

		nx *= 2;
		ny *= 2;
	}

	for( k = 0; k < 4; k++ ) {
		delete[] pts[k];
	}
	delete[] pts;

	return 1;
}
