#include <cmath>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Basis.h"

Basis::Basis( Polygon* _poly, int _order, double* _origin, double _dx, double _dy ) {
	int        i, j, k, xPower, yPower;
	double     fac1, fac2, weight;
	Triangle*  tri;

	poly = _poly;

	order = _order;

	origin[0] = _origin[0];
	origin[1] = _origin[1];

	nFuncs = order*order;

	dxInv = 1.0/_dx;
	dyInv = 1.0/_dy;
	aInv = 1.0/poly->Area();

	ci = new double[nFuncs];
	for( i = 0; i < nFuncs; i++ ) {
		ci[i] = 0.0;
	}
	mean = new double[nFuncs];
	scale = new double[nFuncs];

	/* initialize the mean and scale components */
	scale[0] = 1.0;
	mean[0] = 0.0;
	for( i = 1; i < nFuncs; i++ ) {
		xPower = i%order;
		yPower = i/order;
		fac1 = fac2 = 1.0;

		for( j = 1; j < xPower; j++ ) {
			fac1 *= j;
		}
		for( k = 1; k < yPower; k++ ) {
			fac2 *= k;
		}

		/* normalise coefficients to improve condition number of the matrix */
		scale[i] = pow( dxInv, xPower )*pow( dyInv, yPower )/fac1/fac2;
		
		/* remove mean component so higher order basis functions are massless */
		for( j = 0; j < poly->n; j++ ) {
			for( k = 0; k < poly->tris[j]->nQuadPts; k++ ) {
				tri = poly->tris[j];
				weight = tri->wi[k]*tri->Area();
				mean[i] += weight*pow( tri->qi[k][0] - origin[0], xPower )*pow( tri->qi[k][1] - origin[1], yPower );
			}
		}
		mean[i] *= aInv;
	}
}

Basis::~Basis() {
	delete[] ci;
	delete[] mean;
	delete[] scale;
}

double Basis::EvalIJ( double* pt, int i ) {
	int			xPower = i%order;
	int			yPower = i/order;

	if( i == 0 ) {
		return 1.0;
	}

	return scale[i]*( pow( pt[0] - origin[0], xPower )*pow( pt[1] - origin[1], yPower ) - mean[i] );
}

double Basis::EvalDerivIJ( double* pt, int i, int dim ) {
	int			xPower = i%order;
	int			yPower = i/order;	
	int			coeff = ( dim == 0 ) ? xPower : yPower;

	if( i == 0 ) {
		return 0.0;
	}

	if( dim == 0 && xPower > 0 ) {
		xPower--;
	}
	else if( dim == 1 && yPower > 0 ) {
		yPower--;
	}

	return scale[i]*( pow( pt[0] - origin[0], xPower )*pow( pt[1] - origin[1], yPower ) - mean[i] );
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

double Basis::EvalDerivFull( double* pt, int dim ) {
	int 	i;
	double	result	= 0.0;

	for( i = 1; i < nFuncs; i++ ) {
		result += ci[i]*EvalDerivIJ( pt, i, dim );
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

bool Basis::TestMean( double* volErr ) {
	int 		i, j, k;
	double		weight;
	Triangle*	tri;
	bool		ok;

	*volErr = 0.0;

	for( j = 0; j < poly->n; j++ ) {
		tri = poly->tris[j];
		for( k = 0; k < tri->nQuadPts; k++ ) {
			weight = tri->wi[k]*tri->Area()/poly->Area();
			for( i = 1; i < nFuncs; i++ ) {
				*volErr += weight*ci[i]*EvalIJ( tri->qi[k], i );
			}
		}
	}

	ok = fabs( *volErr ) > 1.0e-12 ? false : true;

	return ok;
}
