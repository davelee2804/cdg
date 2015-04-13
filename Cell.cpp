#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"

Cell::Cell( double** _verts, int _n, int _quadOrder, int _nc, double** _coords ) : Polygon( _verts, _n, _quadOrder ) {
	int i;

	nc = _nc;
	coords = new double*[nc];

	for( i = 0; i < nc; i++ ) {
		coords[i] = new double[2];
		coords[i][0] = _coords[i][0];
		coords[i][1] = _coords[i][1];
	}
}

Cell::~Cell() {
	int i;

	for( i = 0; i < nc; i++ ) {
		delete[] coords[i];
	}
	delete[] coords;
}
