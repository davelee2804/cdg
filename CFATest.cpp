#include <cstdlib>
#include <iostream>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Grid.h"
#include "Basis.h"
#include "Field.h"
#include "CFA.h"

using namespace std;

#define NX 32
#define NY 32

#define QUAD_ORDER 2

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
	double r2 = (p[0] - 0.0)*(p[0] - 0.0) + (p[1] - 0.5)*(p[1] - 0.5);
	if( sqrt( r2 ) < 0.64 ) return exp(-40.0*r2);
	return 0.0;
}

double p1( double* p ) {
	double r2 = (p[0] + 0.5)*(p[0] + 0.5) + (p[1] - 0.0)*(p[1] - 0.0);
	if( sqrt( r2 ) < 0.64 ) return exp(-40.0*r2);
	return 0.0;
}

int main() {
	Grid*	pgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, 1, true );
	Grid*	vgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, 2, false );
	Field*	velx	= new Field( vgrid );
	Field*	vely	= new Field( vgrid );
	Field*	phi		= new Field( pgrid );
	CFA*	cfa		= new CFA( phi, velx, vely );
	int		i, j;
	int		nsteps	= 64*4;
	int		dump	= 1;
	double	dt		= 0.5*M_PI/nsteps;
	Field*	ans		= new Field( pgrid );
	double	err		= 0.0;
	double	norm	= 0.0;
	double	coord[2];

	for( i = 0; i < vgrid->nPolys; i++ ) {
		for( j = 0; j < 4; j++ ) {
			coord[0] = vgrid->minx + (i%vgrid->nx)*vgrid->dx + (j%2)*vgrid->dx;
			coord[1] = vgrid->miny + (i/vgrid->nx)*vgrid->dy + (j/2)*vgrid->dy;
			velx->basis[i]->ci[j] = ux( coord );
			vely->basis[i]->ci[j] = uy( coord );
		}
	}
	for( i = 0; i < pgrid->nPolys; i++ ) {
		coord[0] = pgrid->minx + (i%pgrid->nx)*pgrid->dx + 0.5*vgrid->dx;
		coord[1] = pgrid->miny + (i/pgrid->nx)*pgrid->dy + 0.5*vgrid->dy;
		phi->basis[i]->ci[0] = p0( coord );
		ans->basis[i]->ci[0] = p1( coord );
	}

	pgrid->Write( "pgrid", 2 );
	vgrid->Write( "vgrid", 2 );
	velx->Write( "velx", 0, 2 );
	vely->Write( "vely", 0, 2 );
	phi->Write( "phi", 0, 2 );
	ans->Write( "ans", 0, 2 );

	for( i = 1; i <= nsteps; i++ ) {
		cout << "time step: " << i << "\tdt: " << dt << endl;
		cfa->Advect( dt );
		if( i%dump == 0 ) {
			phi->Write( "phi", i, 2 );
		}
	}

	for( i = 0; i < pgrid->nPolys; i++ ) {
		if( i%NX == 0 || i%NX == NX-1 || i/NX == 0 || i/NX == NY-1 ) {
			continue;
		}
		err += fabs(phi->basis[i]->ci[0] - ans->basis[i]->ci[0]);
		norm += ans->basis[i]->ci[0];
	}
	cout << "error: " << err/norm << endl;
	cout << "analytic mass: " << ans->IntegrateConstant() << endl;
	cout << "numeric mass:  " << phi->IntegrateConstant() << endl;
	cout << "mass loss:     " << 1.0 - phi->IntegrateConstant()/ans->IntegrateConstant() << endl;

	delete pgrid;
	delete vgrid;
	delete velx;
	delete vely;
	delete phi;
	delete cfa;

	return 1;
}
