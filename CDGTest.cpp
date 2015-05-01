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
#include "LinAlg.h"
#include "CDG.h"

using namespace std;

#define NX 32
#define NY 32

#define QUAD_ORDER 3
#define BASIS_ORDER 3

double ux( double* p ) {
	double radius	= sqrt( p[0]*p[0] + p[1]*p[1] );
	double theta 	= atan2( p[1], p[0] );
	double result	= -radius*sin( theta );
	return result;
}

double uy( double* p ) {
	double radius	= sqrt( p[0]*p[0] + p[1]*p[1] );
	double theta 	= atan2( p[1], p[0] );
	double result	= +radius*cos( theta );
	return result;
}

double p0( double* p ) {
	double r2 = (p[0] + 0.0)*(p[0] + 0.0) + (p[1] - 0.5)*(p[1] - 0.5);
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

double p1( double* p ) {
	double r2 = (p[0] + 0.5)*(p[0] + 0.5) + (p[1] - 0.0)*(p[1] - 0.0);
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

int main() {
	Grid*		vgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, 2, false );
	Grid*		pgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
	Field*		velx	= new Field( vgrid );
	Field*		vely	= new Field( vgrid );
	Field*		phi		= new Field( pgrid );
	CDG*		cdg;
	int			i, j;
	int			nsteps	= 64*4;
	int			dump	= 1;
	double		dt		= 0.5*M_PI/nsteps;
	Field*		ans		= new Field( pgrid );
	double		err		= 0.0;
	double		norm	= 0.0;
	double		val_n, val_a;

	for( i = 0; i < vgrid->nCells; i++ ) {
		for( j = 0; j < 4; j++ ) {
			velx->basis[i]->ci[j] = ux( vgrid->cells[i]->coords[j] );
			vely->basis[i]->ci[j] = uy( vgrid->cells[i]->coords[j] );
		}
	}

	/* set up the final solution */
	cdg	= new CDG( ans, velx, vely );
	cdg->InitBetaIJInv( p1 );
	delete cdg;

	/* set up the actual solver */
	cdg = new CDG( phi, velx, vely );
	cdg->InitBetaIJInv( p0 );

	pgrid->Write( "pgrid" );
	vgrid->Write( "vgrid" );
	velx->Write( "velx", 0 );
	vely->Write( "vely", 0 );
	ans->Write( "ans", 0 );
	phi->Write( "phi", 0 );

	cout << "volume: " << phi->Integrate() << endl;

	for( i = 1; i <= nsteps; i++ ) {
		cout << "time step: " << i;
		cdg->Advect( dt );
		cout << "\t...done, volume: " << phi->Integrate() << endl;
		if( i%dump == 0 ) {
			phi->Write( "phi", i );
		}
	}

	for( i = 0; i < pgrid->nCells; i++ ) {
		if( i%NX == 0 || i%NX == NX-1 || i/NX == 0 || i/NX == NY-1 ) {
			continue;
		}

		for( j = 0; j < BASIS_ORDER*BASIS_ORDER; j++ ) {
			val_n = phi->basis[i]->EvalFull( pgrid->cells[i]->coords[j] );
			val_a = ans->basis[i]->EvalFull( pgrid->cells[i]->coords[j] );
			err += fabs(val_a - val_n);
			norm += ans->basis[i]->ci[j];
		}
	}
	cout << "error:      " << err/norm << endl;
	cout << "mass loss:  " << 1.0 - phi->IntegrateConstant()/ans->IntegrateConstant() << endl;

	delete cdg;
	delete phi;
	delete velx;
	delete vely;
	delete pgrid;
	delete vgrid;

	return 1;
}
