#include <cstdlib>
#include <iostream>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"
#include "CFA.h"
#include "CDG.h"

using namespace std;

#define NX 8
#define NY 8

#define QUAD_ORDER 2
#define BASIS_ORDER 2

int main() {
	Grid*	pgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
	Grid*	vgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, 2, false );
	Field*	velx	= new Field( vgrid );
	Field*	vely	= new Field( vgrid );
	Field*	phi		= new Field( pgrid );
	CDG*	cdg		= new CDG( phi, velx, vely );
	int		i, j;
	double	dt		= 0.25*2.0/NX/1.0;

	for( i = 0; i < vgrid->nCells; i++ ) {
		for( j = 0; j < vgrid->cells[j]->nc; j++ ) {
			velx->basis[i]->ci[j] = +1.0;
			vely->basis[i]->ci[j] = +1.0;
		}
	}
	for( j = 0; j < pgrid->cells[0]->nc; j++ ) {
		phi->basis[1*NX+1]->ci[j] = 1.0;
	}

	vgrid->Write( "vgrid" );
	pgrid->Write( "pgrid" );
	velx->Write( "velx", 0 );
	vely->Write( "vely", 0 );
	phi->Write( "phi", 0 );

	for( i = 1; i <= 12; i++ ) {
		cout << "time step: " << i << endl;
		cdg->Advect( dt );
		phi->Write( "phi", i );
	}

	delete cdg;
	delete phi;
	delete velx;
	delete vely;
	delete pgrid;
	delete vgrid;

	return 1;
}
