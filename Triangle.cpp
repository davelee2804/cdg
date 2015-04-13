#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"

using namespace std;

Triangle::Triangle( double* x1, double* x2, double* x3, int _order ) {
	a[0] = x1[0];
	a[1] = x1[1];
	b[0] = x2[0];
	b[1] = x2[1];
	c[0] = x3[0];
	c[1] = x3[1];

	order = _order;

	if( order == 1 ) {
		double tmp[1][2];

		nQuadPts = 1;
		wi = new double[nQuadPts];
		wi[0] = 1.0;
		tmp[0][0] = 1.0/3.0;
		tmp[0][1] = 1.0/3.0;
		qi = new double*[nQuadPts];
		qi[0] = new double[2];
		qi[0][0] = a[0]*tmp[0][0] + b[0]*tmp[0][1] + c[0]*(1.0 - tmp[0][0] - tmp[0][1]);
		qi[0][1] = a[1]*tmp[0][0] + b[1]*tmp[0][1] + c[1]*(1.0 - tmp[0][0] - tmp[0][1]);
	}
	else if( order == 2 ) {
		int i;
		double tmp[3][2];

		nQuadPts = 3;

		wi = new double[nQuadPts];
		wi[0] = 1.0/3;
		wi[1] = 1.0/3;
		wi[2] = 1.0/3;

		tmp[0][0] = 1.0/6.0;
		tmp[0][1] = 2.0/3.0;

		tmp[1][0] = 2.0/3.0;
		tmp[1][1] = 1.0/6.0;

		tmp[2][0] = 1.0/6.0;
		tmp[2][1] = 1.0/6.0;

		qi = new double*[nQuadPts];
		for( i = 0; i < nQuadPts; i++ ) {
			/* generate the coordinates of the quadrature points from the barycentric coordinates */
			qi[i] = new double[2];
			qi[i][0] = a[0]*tmp[i][0] + b[0]*tmp[i][1] + c[0]*(1.0 - tmp[i][0] - tmp[i][1]);
			qi[i][1] = a[1]*tmp[i][0] + b[1]*tmp[i][1] + c[1]*(1.0 - tmp[i][0] - tmp[i][1]);
		}
	}
	else {
		cerr << "triangle quadrature order: " << order << "not implemented" << endl;
		abort();
	}

	ab = new Edge( a, b );
	bc = new Edge( b, c );
	ca = new Edge( c, a );
}

Triangle::~Triangle() {
	int i;

	for( i = 0; i < nQuadPts; i++ ) {
		delete[] qi[i];
	}
	delete[] qi;
	delete[] wi;

	delete ab;
	delete bc;
	delete ca;
}

double Triangle::Area() {
	return 0.5*fabs( (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]) );
}

double Triangle::CrossProduct( int origin ) {
	double 	p1[2], p2[2], p3[2], x12, y12, x13, y13;

	if( origin == 1 ) {
		p1[0] = a[0]; p2[0] = b[0]; p3[0] = c[0];
		p1[1] = a[1]; p2[1] = b[1]; p3[1] = c[1];
	}
	else if( origin == 2 ) {
		p1[0] = b[0]; p2[0] = c[0]; p3[0] = a[0];
		p1[1] = b[1]; p2[1] = c[1]; p3[1] = a[1];
	}
	else if( origin == 3 ) {
		p1[0] = c[0]; p2[0] = a[0]; p3[0] = b[0];
		p1[1] = c[1]; p2[1] = a[1]; p3[1] = b[1];
	}
	else {
		cerr << "invalid origin vertex index (must be 1-3)" << endl;
		abort();
	}

	x12 = p2[0] - p1[0];
	y12 = p2[1] - p1[1];
	x13 = p3[0] - p1[0];
	y13 = p3[1] - p1[1];

	return x12*y13 - y12*x13;
}