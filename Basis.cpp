#include <cmath>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Basis.h"

//#define BASIS_TEST 1

Basis::Basis( int _order, double* _origin, Polygon* _poly ) {
	int i;

	order = _order;

	origin[0] = _origin[0];
	origin[1] = _origin[1];

	poly = _poly;

	ci = new double[order*order];
	for( i = 0; i < order*order; i++ ) {
		ci[i] = 0.0;
	}
}

Basis::~Basis() {
	delete[] ci;
}

double Basis::EvalIJ( double* pt, int i ) {
#ifdef BASIS_TEST
	if( !poly->IsInside( pt ) ) {
		cerr << "ERROR: basis function to be evaluated at point outside corresponding polygon..." << endl;
		abort(); 
	}
#endif

	return pow( pt[0] - origin[0], i%order )*pow( pt[1] - origin[1], i/order );
}

double Basis::EvalConst( double* pt ) {
	int 	i;
	double	result	= 0.0;

	for( i = 0; i < order*order; i++ ) {
		result += EvalIJ( pt, i );
	}
	return result;
}

double Basis::EvalFull( double* pt ) {
	int 	i;
	double	result	= 0.0;

	for( i = 0; i < order*order; i++ ) {
		result += ci[i]*EvalIJ( pt, i );
	}
	return result;
}
