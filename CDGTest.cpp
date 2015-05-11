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
#include "LinAlg.h"
#include "CDG.h"
#include "Limiter.h"

//#define DISTORT_MESH

using namespace std;

#define NX 32
#define NY 32

#define QUAD_ORDER 3
#define BASIS_ORDER 2

double ux( double* p ) {
	return -p[1];
}

double uy( double* p ) {
	return +p[0];
}

double p0( double* p ) {
	double	xo	= 0.5*cos( 0.25*M_PI );
	double	yo	= 0.5*sin( 0.25*M_PI );
	double	r2	= ( p[0] - xo )*( p[0] - xo ) + ( p[1] - yo )*( p[1] - yo );
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

double p1( double* p ) {
	double	xo	= 0.5*cos( 0.75*M_PI );
	double	yo	= 0.5*sin( 0.75*M_PI );
	double	r2	= ( p[0] - xo )*( p[0] - xo ) + ( p[1] - yo )*( p[1] - yo );
	if( sqrt( r2 ) < 0.40 ) return exp(-40.0*r2);
	return 0.0;
}

void TestQuadArea( Polygon* poly ) {
	double		side[2], p[2], q[2], s2[4], p2, q2, det, a1, a2;
	int			i;

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

int main() {
	Grid*		pgrid 	= new Grid( NX, NY, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER, true );
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
	double		coord[2];
	Limiter*	lim		= new Limiter( phi );

#ifdef DISTORT_MESH
	srand( 7919 );
	for( i = 0; i < pgrid->nVerts; i++ ) {
		if( i%(pgrid->nx + 1) == 0 || i%(pgrid->nx + 1) == pgrid->nx || i/(pgrid->nx + 1) == 0 || i/(pgrid->nx + 1) == pgrid->ny ) {
			continue;
		}

		pgrid->verts[i][0] += 0.1*pgrid->dx*( (2.0*rand())/RAND_MAX - 1.0 );
		pgrid->verts[i][1] += 0.1*pgrid->dy*( (2.0*rand())/RAND_MAX - 1.0 );
	}
	pgrid->UpdateEdges();
	pgrid->UpdatePolys();
	pgrid->UpdateTris();
	phi->UpdateBasis();
	ans->UpdateBasis();
	for( i = 0; i < pgrid->nPolys; i++ ) {
		TestQuadArea( pgrid->polys[i] );
	}
#endif

	/* set up the final solution */
	cdg	= new CDG( ans, NULL, NULL, ux, uy );
	cdg->InitBetaIJInv( p1 );
	delete cdg;

	/* set up the actual solver */
	cdg = new CDG( phi, NULL, NULL, ux, uy );
	cdg->InitBetaIJInv( p0 );

	pgrid->Write( "pgrid", BASIS_ORDER );
	ans->Write( "ans", 0, BASIS_ORDER );
	phi->Write( "phi", 0, BASIS_ORDER );

	cout << "volume: " << phi->Integrate() << endl;

	for( i = 1; i <= nsteps; i++ ) {
		cout << "time step: " << i;
		cdg->Advect( dt );
		//lim->Apply();

		cout << "\t...done, volume: " << phi->Integrate() << endl;
		if( i%dump == 0 ) {
			phi->Write( "phi", i, BASIS_ORDER );
		}
	}

	for( i = 0; i < pgrid->nPolys; i++ ) {
		if( i%NX == 0 || i%NX == NX-1 || i/NX == 0 || i/NX == NY-1 ) {
			continue;
		}

		for( j = 0; j < BASIS_ORDER*BASIS_ORDER; j++ ) {
			coord[0] = pgrid->minx + (i%pgrid->nx)*pgrid->dx + (0.5 + j%BASIS_ORDER)/BASIS_ORDER*pgrid->dx;
			coord[1] = pgrid->miny + (i/pgrid->nx)*pgrid->dy + (0.5 + j/BASIS_ORDER)/BASIS_ORDER*pgrid->dy;
			val_n = phi->basis[i]->EvalFull( coord );
			val_a = ans->basis[i]->EvalFull( coord );
			err += fabs(val_a - val_n);
			norm += ans->basis[i]->ci[j];
		}
	}
	cout << "L_1 error:  " << err/norm << endl;
	cout << "L_2 error:  " << phi->L2Error( ans ) << endl;
	cout << "mass loss:  " << 1.0 - phi->Integrate()/ans->Integrate() << endl;

	delete lim;
	delete cdg;
	delete phi;
	delete pgrid;

	return 1;
}
