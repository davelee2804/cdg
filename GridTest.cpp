#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Grid.h"

using namespace std;

#define QUAD_ORDER 2
#define BASIS_ORDER 2

double Rand( double min, double max ) { return ((double)(max - min)*rand())/RAND_MAX + min; }

bool PointDiff( double* p1, double* p2, double tol ) {
	if( sqrt( (p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p1[1])*(p1[1] - p1[1]) ) > tol ) {
		return false;
	}
	return true;
}

int main() {
	int			nx			= 6;
	int			ny			= 4;
	int 		minx		= -1.0;
	int			miny		= -1.0;
	int			maxx		= +1.0;
	int			maxy		= +1.0;
	Grid*		grid		= new Grid( nx, ny, minx, miny, maxx, maxy, QUAD_ORDER, BASIS_ORDER, false );
	double		point[2]	= { 0.0, 0.0 };
	Polygon*	testPoly	= NULL;
	Polygon*	testPoly2	= NULL;
	Polygon*	testPoly3	= NULL;
	int			i, j;
	bool		hasEdge;

	for( i = 0; i < 100; i++ ) {
		point[0] = Rand( minx, maxx );
		point[1] = Rand( miny, maxy );

		//point->Print();

		testPoly = grid->FindPoly( point );
	}

	for( i = 0; i < grid->nPolys; i++ ) {
		testPoly = grid->polys[i];
		if( !PointDiff( testPoly->verts[0], testPoly->edges[0]->v1, 0.001 ) ) {
			cout << "ERROR [0]" << endl;
		}
		if( !PointDiff( testPoly->verts[0], testPoly->edges[2]->v1, 0.001 ) ) {
			cout << "ERROR [1]" << endl;
		}
		if( !PointDiff( testPoly->verts[1], testPoly->edges[0]->v2, 0.001 ) ) {
			cout << "ERROR [2]" << endl;
		}
		if( !PointDiff( testPoly->verts[1], testPoly->edges[3]->v1, 0.001 ) ) {
			cout << "ERROR [3]" << endl;
		}
		if( !PointDiff( testPoly->verts[2], testPoly->edges[1]->v2, 0.001 ) ) {
			cout << "ERROR [4]" << endl;
		}
		if( !PointDiff( testPoly->verts[2], testPoly->edges[3]->v2, 0.001 ) ) {
			cout << "ERROR [5]" << endl;
		}
		if( !PointDiff( testPoly->verts[3], testPoly->edges[1]->v1, 0.001 ) ) {
			cout << "ERROR [6]" << endl;
		}
		if( !PointDiff( testPoly->verts[3], testPoly->edges[2]->v2, 0.00 ) ) {
			cout << "ERROR [7]" << endl;
		}
	}

	for( i = 0; i < grid->nEdges; i++ ) {
		testPoly2 = testPoly3 = NULL;
		grid->GetEdgePolys( i, testPoly2, testPoly3 );
		if( testPoly2 && testPoly3 ) {
			hasEdge = false;
			for( j = 0; j < 4; j++ )
				if( testPoly2->edges[j] == grid->edges[i] )
					hasEdge = true;

			if( !hasEdge )
				cout << "ERROR! [8]" << endl;

			hasEdge = false;
			for( j = 0; j < 4; j++ )
				if( testPoly3->edges[j] == grid->edges[i] )
					hasEdge = true;

			if( !hasEdge )
				cout << "ERROR! [9]" << endl;
		}
	}

	delete grid;

	return 1;
}
