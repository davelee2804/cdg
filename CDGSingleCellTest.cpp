#include <cstdlib>
#include <iostream>
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

#define NX 5
#define NY 5

#define QUAD_ORDER 3
#define BASIS_ORDER 2

double p0( double* x ) {
	if( x[0] > -0.6001 && x[0] < +0.6001 && x[1] > -0.6001 && x[1] < +0.6001 ) return 0.25*x[0];

	return 0.0;
}

double ux( double* x ) { return -0.1; }
double uy( double* x ) { return -0.1; }

int main() {
	Grid*	grid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
	Field*	phi		= new Field( grid );
	CDG*	cdg		= new CDG( phi, NULL, NULL, ux, uy );
	int		i;
	double	dt		= 0.1*(4.0/NX)/0.1;

	cdg->InitBetaIJInv( p0 );

	grid->WriteTris( "pgrid" );
	phi->WriteBasis( "phi", 0 );
	phi->WriteDeriv( "phi", 0, 0 );
	phi->WriteDeriv( "phi", 0, 1 );

	for( i = 1; i <= 1; i++ ) {
		cout << "time step: " << i << endl;
		cdg->Advect( dt );
		phi->WriteBasis( "phi", i );
		phi->WriteDeriv( "phi", i, 0 );
		phi->WriteDeriv( "phi", i, 1 );
	}

	delete cdg;
	delete phi;
	delete grid;

	return 1;
}
