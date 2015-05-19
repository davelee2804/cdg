#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"
#include "LinAlg.h"
#include "CFA.h"
#include "CDG.h"

using namespace std;

#define QUAD_ORDER 6
#define BASIS_ORDER 3

double p1( double* p ) {
	double	xo	= 0.5*cos( 0.75*M_PI );
	double	yo	= 0.5*sin( 0.75*M_PI );
	double	r2	= ( p[0] - xo )*( p[0] - xo ) + ( p[1] - yo )*( p[1] - yo );
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

void WriteDifference( Field* phi, Field* ans, int nx ) {
	Grid*		grid		= phi->grid;
	char 		filename[50];
	ofstream 	file;
	int 		i, j, k, l;
	int 		n 			= 2;
	double		point[2], val;

	sprintf( filename, "output/diff_%.3u.txt", nx );
	file.open( filename );
	for( i = 0; i < grid->ny; i++ ) {
		for( j = 0; j < n; j++ ) {
			for( k = 0; k < grid->nx; k++ ) {
				for( l = 0; l < n; l++ ) {
					point[0] = grid->minx + k*grid->dx + (0.5 + l)*grid->dx/n;
					point[1] = grid->miny + i*grid->dy + (0.5 + j)*grid->dy/n;
					val = ans->basis[i*grid->nx+k]->EvalFull( point ) - phi->basis[i*grid->nx+k]->EvalFull( point );
					file << val << endl;
				}
			}
		}
		file << endl;
	}
	file.close();
}

int main( int argc, char** argv ) {
	int			nx, i;
	Grid*		grid;
	Field*		phi;
	char		filename[80];
	double		l1[3], l2[3], l1_norm[3], l2_norm[3];

	i = 0;
	for( nx = 64; nx < 512; nx *= 2 ) {
		cout << "nx: " << nx << endl;

		grid = new Grid( nx, nx, -2.0, -2.0, +2.0, +2.0, QUAD_ORDER, BASIS_ORDER, true );
		phi	 = new Field( grid );

		sprintf( filename, "input/Q%uB%u_basis_%.3u.txt", QUAD_ORDER, BASIS_ORDER, nx );
		phi->ReadBasis( filename, 0 );

		l1_norm[i] = phi->L1Error( p1, true );
		l2_norm[i] = phi->L2Error( p1, true );
		l1[i] = phi->L1Error( p1, false );
		l2[i] = phi->L2Error( p1, false );

		i++;

		delete phi;
		delete grid;
	}

	cout << "\nL_1 (un-normalized): ";
	for( i = 0; i < 3; i++ ) { cout << l1[i] << "\t"; }
	cout << endl;

	cout << "\nL_2 (un-normalized): ";
	for( i = 0; i < 3; i++ ) { cout << l2[i] << "\t"; }
	cout << endl;

	cout << "\nL_1 (normalized):    ";
	for( i = 0; i < 3; i++ ) { cout << l1_norm[i] << "\t"; }
	cout << endl;

	cout << "\nL_2 (normalized):    ";
	for( i = 0; i < 3; i++ ) { cout << l2_norm[i] << "\t"; }
	cout << endl;

	return 1;
}
