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

using namespace std;

#define QUAD_ORDER 3
#define BASIS_ORDER 2

double p1( double* p ) {
	double	xo	= 0.5*cos( 0.75*M_PI );
	double	yo	= 0.5*sin( 0.75*M_PI );
	double	r2	= ( p[0] - xo )*( p[0] - xo ) + ( p[1] - yo )*( p[1] - yo );
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

double ErrorPerPoly( Field* phi, Func* f ) {
	Grid*		grid	= phi->grid;
	Polygon*	poly;
	double		error 	= 0.0;
	int			i;

	for( i = 0; i < grid->nPolys; i++ ) {
		poly = grid->polys[i];
		error += grid->dx*grid->dy*fabs( phi->basis[i]->EvalFull( poly->origin ) - f( poly->origin ) );
	}
	return error;
}

int main( int argc, char** argv ) {
	int			nx, i;
	Grid*		grid;
	Field*		phi;
	char		filename[80];
	double		l1[3], l2[3], l1_norm[3], l2_norm[3], epp[3];

	i = 0;
	for( nx = 32; nx < 256; nx *= 2 ) {
		cout << "nx: " << nx << endl;

		grid 	= new Grid( nx, nx, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
		phi		= new Field( grid );

		sprintf( filename, "input/Q%uB%u_basis_%.3u.txt", QUAD_ORDER, BASIS_ORDER, nx );
		phi->ReadBasis( filename );

		sprintf( filename, "phi_%.3u", nx );
		phi->Write( filename, 256, 1 );

		sprintf( filename, "pgrid_%.3u", nx );
		grid->Write( filename, 1 );

		l1_norm[i] = phi->L1Error( p1, true );
		l2_norm[i] = phi->L2Error( p1, true );
		l1[i] = phi->L1Error( p1, false );
		l2[i] = phi->L2Error( p1, false );
		epp[i] = ErrorPerPoly( phi, p1 );

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

	cout << "\nerror per cell:      ";
	for( i = 0; i < 3; i++ ) { cout << epp[i] << "\t"; }
	cout << endl;

	return 1;
}
