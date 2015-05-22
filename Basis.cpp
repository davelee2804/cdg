#include <cmath>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Basis.h"

Basis::Basis( Polygon* _poly, int _order, double* _origin, double _dx, double _dy ) {
	int i;

	poly = _poly;

	order = _order;

	origin[0] = _origin[0];
	origin[1] = _origin[1];

	if( order == 0 ) {
		nFuncs = 1;
	}
	else if( order == 1 ) {
		nFuncs = 3;
	}
	else if( order == 2 ) {
		nFuncs = 6;
	}
	else {
		cout << "ERROR: basis order " << order << " not yet implemented..." << endl;
		abort(); 
	}
	xPower = new int[nFuncs];
	yPower = new int[nFuncs];

	xPower[0] = yPower[0] = 0;
	if( order == 1 ) {
		xPower[1] = 1;
		xPower[2] = 0;
		yPower[1] = 0;
		yPower[2] = 1;
	}
	if( order == 2 ) {
		xPower[1] = 1;
		xPower[2] = 0;
		xPower[3] = 2;
		xPower[4] = 1;
		xPower[5] = 0;

		yPower[1] = 0;
		yPower[2] = 1;
		yPower[3] = 0;
		yPower[4] = 1;
		yPower[5] = 2;
	}

	dxInv = 1.0/_dx;
	dyInv = 1.0/_dy;

	ci = new double[nFuncs];
	for( i = 0; i < nFuncs; i++ ) {
		ci[i] = 0.0;
	}
	scale = new double[nFuncs];
	mean = new double[nFuncs];

	Init();
}

Basis::~Basis() {
	delete[] ci;
	delete[] scale;
	delete[] mean;
	delete[] xPower;
	delete[] yPower;
}

void Basis::Init() {
	int        i, j, k;
	double     fac1, fac2, weight;
	Triangle*  tri;

	/* initialize the mean and scale components */
	scale[0] = 1.0;
	mean[0] = 0.0;
	for( i = 1; i < nFuncs; i++ ) {
		fac1 = fac2 = 1.0;

		for( j = 1; j <= xPower[i]; j++ ) {
			fac1 *= j;
		}
		for( k = 1; k <= yPower[i]; k++ ) {
			fac2 *= k;
		}

		/* normalise coefficients to improve condition number of the matrix */
		scale[i] = pow( dxInv, xPower[i] )*pow( dyInv, yPower[i] )/fac1/fac2;
		
		/* remove mean component so higher order basis functions are massless */
		mean[i] = 0.0;
		for( j = 0; j < poly->n; j++ ) {
			tri = poly->tris[j];
			for( k = 0; k < tri->nq; k++ ) {
				weight = tri->wi[k]*tri->area;
				mean[i] += weight*pow( tri->qi[k][0] - origin[0], xPower[i] )*pow( tri->qi[k][1] - origin[1], yPower[i] );
			}
		}
		mean[i] /= poly->Area();
	}
}

double Basis::EvalIJ( double* pt, int i ) {
	if( i == 0 ) {
		return 1.0;
	}

	return scale[i]*( pow( pt[0] - origin[0], xPower[i] )*pow( pt[1] - origin[1], yPower[i] ) - mean[i] );
}

double Basis::EvalDerivIJ( double* pt, int i, int dim ) {
	int xp    = xPower[i];
	int yp    = yPower[i];
	int coeff = ( dim == 0 ) ? xp : yp;

	if( i == 0 ) {
		return 0.0;
	}

	if( dim == 0 && xp > 0 ) {
		xp--;
	}
	else if( dim == 1 && yp > 0 ) {
		yp--;
	}

	return scale[i]*( coeff*pow( pt[0] - origin[0], xp )*pow( pt[1] - origin[1], yp ) - mean[i] );
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
	int        i, j, k;
	double     weight;
	Triangle*  tri;
	bool       ok;

	*volErr = 0.0;

	for( j = 0; j < poly->n; j++ ) {
		tri = poly->tris[j];
		for( k = 0; k < tri->nq; k++ ) {
			weight = tri->wi[k]*tri->area;
			for( i = 1; i < nFuncs; i++ ) {
				*volErr += weight*ci[i]*EvalIJ( tri->qi[k], i );
			}
		}
	}

	ok = fabs( *volErr ) > 1.0e-12 ? false : true;

	return ok;
}
