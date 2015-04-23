#include <cmath>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Cell.h"
#include "Basis.h"

Basis::Basis( int _order, double* _origin ) {
	int i;

	order = _order;

	origin[0] = _origin[0];
	origin[1] = _origin[1];

	nFuncs = order*order;

	ci = new double[nFuncs];
	for( i = 0; i < nFuncs; i++ ) {
		ci[i] = 0.0;
	}
}

Basis::~Basis() {
	delete[] ci;
}

double Basis::EvalIJ( double* pt, int i ) {
	return pow( pt[0] - origin[0], i%order )*pow( pt[1] - origin[1], i/order );
}

double Basis::EvalConst( double* pt ) {
	int 	i;
	double	result	= 0.0;

	for( i = 0; i < nFuncs; i++ ) {
		result += EvalIJ( pt, i );
	}
	return result;
}

double Basis::EvalFull( double* pt ) {
	int 	i;
	double	result	= 0.0;

	for( i = 0; i < nFuncs; i++ ) {
		result += ci[i]*EvalIJ( pt, i );
	}
	return result;
}

double Basis::EvalWithCoeffs( double* pt, double* coeffs ) {
	int 	i;
	double	result	= 0.0;

	for( i = 0; i < nFuncs; i++ ) {
		result += coeffs[i]*EvalIJ( pt, i );
	}
	return result;
}
