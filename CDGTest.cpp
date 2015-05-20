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
#include "CFA.h"
#include "LinAlg.h"
#include "CDG.h"
#include "Limiter.h"

//#define DISTORT_MESH

using namespace std;

#define NX 64
#define NY 64

#define QUAD_ORDER 3
#define BASIS_ORDER 2

double ux( double* p ) {
	return -p[1];
}

double uy( double* p ) {
	return +p[0];
}

double p0( double* p ) {
	double xo = 0.5*cos( 0.25*M_PI );
	double yo = 0.5*sin( 0.25*M_PI );
	double r2 = ( p[0] - xo )*( p[0] - xo ) + ( p[1] - yo )*( p[1] - yo );
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

double p1( double* p ) {
	double xo = 0.5*cos( 0.75*M_PI );
	double yo = 0.5*sin( 0.75*M_PI );
	double r2 = ( p[0] - xo )*( p[0] - xo ) + ( p[1] - yo )*( p[1] - yo );
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

void TestQuadArea( Polygon* poly ) {
	double side[2], p[2], q[2], s2[4], p2, q2, det, a1, a2;
	int    i;

	for( i = 0; i < 4; i++ ) {
		side[0] = poly->verts[(i+1)%4][0] - poly->verts[i][0];
		side[1] = poly->verts[(i+1)%4][1] - poly->verts[i][1];
		s2[i] = side[0]*side[0] + side[1]*side[1];
	}
	det = s2[1] + s2[3] - s2[0] - s2[2];

	p[0] = poly->verts[2][0] - poly->verts[0][0];
	p[1] = poly->verts[2][1] - poly->verts[0][1];
	q[0] = poly->verts[3][0] - poly->verts[1][0];
	q[1] = poly->verts[3][1] - poly->verts[1][1];
	p2 = p[0]*p[0] + p[1]*p[1];
	q2 = q[0]*q[0] + q[1]*q[1];

	a1 = 0.25*sqrt( 4.0*p2*q2 - det*det );
	a2 = poly->Area();

	if( fabs( a1 - a2 ) > 1.0e-8 ) {
		cout << "ERROR: area doesn't match expected... " << fabs( a1 - a2 ) << endl;
	}
}

void WriteVelocity( Grid* grid ) {
	ofstream    file;
	char        filename[80];
	int         i;

	sprintf( filename, "output/vgrid.x.txt" );
	file.open( filename );
	for( i = 0; i <= grid->nx; i++ ) {
		file << grid->verts[i][0] << endl;
	}
	file.close();

	sprintf( filename, "output/vgrid.y.txt" );
	file.open( filename );
	for( i = 0; i <= grid->ny; i++ ) {
		file << grid->verts[i*(grid->nx+1)][1] << endl;
	}
	file.close();

	sprintf( filename, "output/velx.0000.txt" );
	file.open( filename );
	for( i = 0; i < grid->nVerts; i++ ) {
		file << ux( grid->verts[i] ) << endl;
	}
	file.close();

	sprintf( filename, "output/vely.0000.txt" );
	file.open( filename );
	for( i = 0; i < grid->nVerts; i++ ) {
		file << uy( grid->verts[i] ) << endl;
	}
	file.close();
}

int main() {
	Grid*     grid      = new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
	Field*    phi       = new Field( grid );
	CDG*      cdg;
	int       i;
	int       start     = 0;
	int       nsteps    = 64*4;
	int       dump      = 8;
	double    dt        = 0.5*M_PI/nsteps;
	Field*    ans       = new Field( grid );
	Limiter*  lim       = new Limiter( phi );

#ifdef DISTORT_MESH
	srand( 7919 );
	for( i = 0; i < grid->nVerts; i++ ) {
		if( i%(grid->nx + 1) == 0 || i%(grid->nx + 1) == grid->nx || i/(grid->nx + 1) == 0 || i/(grid->nx + 1) == grid->ny ) {
			continue;
		}

		grid->verts[i][0] += 0.1*grid->dx*( (2.0*rand())/RAND_MAX - 1.0 );
		grid->verts[i][1] += 0.1*grid->dy*( (2.0*rand())/RAND_MAX - 1.0 );
	}
	grid->UpdateEdges();
	grid->UpdatePolys();
	grid->UpdateTris();
	phi->UpdateBasis();
	ans->UpdateBasis();
	for( i = 0; i < grid->nPolys; i++ ) {
		TestQuadArea( grid->polys[i] );
	}
#endif

	/* set up the final solution */
	cdg	= new CDG( ans, NULL, NULL, ux, uy );
	cdg->InitBetaIJInv( p1 );
	delete cdg;

	/* set up the actual solver */
	if( !start ) {
		cdg = new CDG( phi, NULL, NULL, ux, uy );
		cdg->InitBetaIJInv( p0 );

		grid->WriteTris( "pgrid" );
		ans->WriteBasis( "ans", 0 );
		phi->WriteBasis( "phi", 0 );
		phi->WriteDeriv( "phi", 0, 0 );
		phi->WriteDeriv( "phi", 0, 1 );
		WriteVelocity( grid );
	}
	else {
		phi->ReadBasis( "phi", start );
		cdg = new CDG( phi, NULL, NULL, ux, uy );
		cdg->InitBetaIJInv( NULL );
	}
	cout << "volume: " << phi->Integrate() << endl;

	for( i = start + 1; i <= nsteps; i++ ) {
		cout << "time step: " << i;
		cdg->Advect( dt );
		//lim->Apply();

		cout << "\t...done, volume: " << phi->Integrate() << endl;
		if( i%dump == 0 ) {
			phi->WriteBasis( "phi", i );
			phi->WriteDeriv( "phi", i, 0 );
			phi->WriteDeriv( "phi", i, 1 );
		}
	}

	cout << "L_1 error:  " << phi->L1Error( p1, false ) << endl;
	cout << "L_2 error:  " << phi->L2Error( p1, false ) << endl;
	cout << "mass loss:  " << 1.0 - phi->Integrate()/ans->Integrate() << endl;

	delete lim;
	delete cdg;
	delete phi;
	delete grid;

	return 1;
}
