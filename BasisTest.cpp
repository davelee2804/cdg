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
#include "LinAlg.h"
#include "CDG.h"

using namespace std;

#define QUAD_ORDER 2
#define BASIS_ORDER 2

double func( double* x ) {
	return cos( 0.5*M_PI*x[0] )*cos( 0.5*M_PI*x[1] );
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
	Cell* 		cell;
	Polygon*	poly;
	Triangle*	tri;
	Basis*		basis;
	double		ans			= 16.0/M_PI/M_PI;
	double		vol;
	CDG*		cdg;
	Field*		velx;
	Field*		vely;
	double*		beta_ij_2;
	double*		cj;
	double**	verts;
	double**	pts;
	double		origin[2]	= { 0.0, 0.0 };
	double		beta_ij[BASIS_ORDER*BASIS_ORDER*BASIS_ORDER*BASIS_ORDER];
	double		betaInv_ij[BASIS_ORDER*BASIS_ORDER*BASIS_ORDER*BASIS_ORDER];
	double		weight, *coord, fj[BASIS_ORDER*BASIS_ORDER];

	verts = new double*[4];
	for( i = 0; i < 4; i++ ) {
		verts[i] = new double[2];
	}
	pts = new double*[BASIS_ORDER*BASIS_ORDER];
	for( i = 0; i < BASIS_ORDER*BASIS_ORDER; i++ ) {
		pts[i] = new double[2];
	}

	verts[0][0] = -1.0;
	verts[0][1] = -1.0;
	verts[1][0] = -1.0;
	verts[1][1] = +1.0;
	verts[2][0] = +1.0;
	verts[2][1] = +1.0;
	verts[3][0] = +1.0;
	verts[3][1] = -1.0;
	for( i = 0; i < 4; i++ ) {
		Distort( verts[i] );
	}

	poly = new Polygon( verts, 4, QUAD_ORDER );
	basis = new Basis( poly, BASIS_ORDER, origin, 2.0, 2.0 );

	for( j = 0; j < BASIS_ORDER*BASIS_ORDER*BASIS_ORDER*BASIS_ORDER; j++ ) {
		beta_ij[j] = 0.0;
		betaInv_ij[j] = 0.0;
	}
	for( k = 0; k < poly->n; k++ ) {
		tri = poly->tris[k];
		for( l = 0; l < tri->nQuadPts; l++ ) {
			for( j = 0; j < BASIS_ORDER*BASIS_ORDER; j++ ) {
				for( i = 0; i < BASIS_ORDER*BASIS_ORDER; i++ ) {
					weight = tri->wi[l]*tri->Area()/poly->Area();
					coord = tri->qi[l];
					beta_ij[j*BASIS_ORDER*BASIS_ORDER+i] += weight*basis->EvalIJ( coord, i )*basis->EvalIJ( coord, j );
				}
			}
		}
	}
	MatInv( beta_ij, betaInv_ij, BASIS_ORDER*BASIS_ORDER );
	for( j = 0; j < BASIS_ORDER*BASIS_ORDER; j++ ) {
		fj[j] = 0.0;
		for( k = 0; k < poly->n; k++ ) {
			tri = poly->tris[k];
			for( l = 0; l < tri->nQuadPts; l++ ) {
				weight = tri->wi[l]*tri->Area()/poly->Area();
				coord = tri->qi[l];
				/* basis initially set as the spatial values at the cell coordinates */
				fj[j] += weight*basis->EvalIJ( coord, j )*func( coord );
			}
		}
	}
	/* set the initial basis coefficients */
	AXEB( betaInv_ij, fj, basis->ci, BASIS_ORDER*BASIS_ORDER );
	cout << "basis coeffiecients: ";
	for( i = 0; i < BASIS_ORDER*BASIS_ORDER; i++ ) {
		cout << basis->ci[i] << " ";
	}
	cout << endl;

	if( !basis->TestMean( &vol ) ) {
		cout << "ERROR: basis function mean not equal to first component..." << vol << endl;
	}

	delete poly;
	for( i = 0; i < 4; i++ ) {
		delete[] verts[i];
	}
	for( i = 0; i < BASIS_ORDER*BASIS_ORDER; i++ ) {
		delete[] pts[i];
	}
	delete[] verts;
	delete[] pts;

	cout << "testing the basis function matrix inverse..." << endl;
	for( i = 0; i < 8; i++ ) {
		grid = new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		field = new Field( grid );
		velx = new Field( grid );
		vely = new Field( grid );

		cdg = new CDG( field, velx, vely );
		cdg->InitBetaIJInv( func );

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
