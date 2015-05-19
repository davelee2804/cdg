#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Basis.h"
#include "Grid.h"
#include "Field.h"
#include "CFA.h"
#include "LinAlg.h"
#include "CDG.h"

using namespace std;

#define QUAD_ORDER 6
#define BASIS_ORDER 3

double func( double* x ) {
	return cos( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
}

double dxFunc( double* x ) {
	return -0.5*M_PI*sin( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
}

double dyFunc( double* x ) {
	return -0.5*M_PI*cos( 0.5*M_PI*x[0] )*sin( 0.5*M_PI*x[1] );
}

double dxdyFunc( double* x ) {
	return 0.25*M_PI*M_PI*sin( 0.5*M_PI*x[0] )*sin( 0.5*M_PI*x[1] );
}

void Distort( double* pt ) {
	double alpha = 1.1;
	double gamma = 0.8;
	double theta = 0.15*M_PI;

	pt[0] = alpha*cos(theta)*pt[0] - gamma*sin(theta)*pt[1];
	pt[1] = alpha*sin(theta)*pt[0] + gamma*cos(theta)*pt[1];
}

/* test the basis function representation. since the grid is cartesian, use a first order basis:
 * 		c_0 + c_1.x + c_2.y + c_3.xy, 
 * evaluated at the quadrature points */
int main() {
	int			nx			= 1;
	int			ny			= 1;
	int			i, j, k, l;
	Grid*		grid;
	Field*		field;
	Polygon*	poly;
	Triangle*	tri;
	double		ans			= 16.0/M_PI/M_PI;
	double		vol;
	CDG*		cdg;
	int			nBasis		= BASIS_ORDER*BASIS_ORDER;
	double		phi_xn, phi_yn, phi_xa, phi_ya;
	double		err_x, err_y, norm_x, norm_y;
	double		weight;

	grid = new Grid( 1, 1, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
	field = new Field( grid );
	cdg = new CDG( field, NULL, NULL, NULL, NULL );
	for( i = 0; i < 4; i++ ) {
		Distort( grid->polys[0]->verts[i] );
	}
	cdg->InitBetaIJInv( func );

	cout << "basis coeffiecients: ";
	for( i = 0; i < nBasis; i++ ) {
		cout << field->basis[0]->ci[i] << " ";
	}
	cout << endl;

	if( !field->basis[0]->TestMean( &vol ) ) {
		cout << "ERROR: basis function mean not equal to first component..." << vol << endl;
	}

	delete cdg;
	delete field;
	delete grid;

	cout << "testing the basis function matrix inverse..." << endl;
	for( i = 0; i < 8; i++ ) {
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );

		cdg = new CDG( field, NULL, NULL, NULL, NULL );
		cdg->InitBetaIJInv( func );

		vol = field->Integrate();

		cout << fabs( 1.0 - vol/ans ) << endl;

		delete cdg;
		delete field;
		delete grid;

		nx *= 2;
		ny *= 2;
	}

	cout << "testing the basis derivative function..." << endl;
	nx = ny = 1;
	for( i = 0; i < 8; i++ ) {
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );

		cdg = new CDG( field, NULL, NULL, NULL, NULL );
		cdg->InitBetaIJInv( func );

		err_x = err_y = norm_x = norm_y = 0.0;

		for( j = 0; j < grid->nPolys; j++ ) {
			poly = grid->polys[j];
			for( k = 0; k < poly->n; k++ ) {
				tri = poly->tris[k];
				for( l = 0; l < tri->nq; l++ ) {
					weight = tri->wi[l]*tri->area;

					phi_xn = field->basis[j]->EvalDerivFull( tri->qi[l], 0 );
					phi_yn = field->basis[j]->EvalDerivFull( tri->qi[l], 1 );
					phi_xa = dxFunc( tri->qi[l] );
					phi_ya = dyFunc( tri->qi[l] );

					err_x += weight*fabs(phi_xa - phi_xn);
					err_y += weight*fabs(phi_ya - phi_yn);
					norm_x += weight*fabs( phi_xa );
					norm_y += weight*fabs( phi_ya );
				}
			}
		}

		cout << err_x/norm_x << "\t" << err_y/norm_y << endl;

		delete cdg;
		delete field;
		delete grid;

		nx *= 2;
		ny *= 2;
	}

	return 1;
}
