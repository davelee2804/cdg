#include <iostream>
#include <cmath>

#include "Edge.h"

using namespace std;

Edge::Edge( double* x1, double* x2 ) {
	v1[0] = x1[0];
	v1[1] = x1[1];
	v2[0] = x2[0];
	v2[1] = x2[1];

	Init();
}

Edge::~Edge() {
}

void Edge::Init() {
	dx = v1[0] - v2[0];
	dy = v1[1] - v2[1];

	x0 = 0.5*(v1[0] + v2[0]);
	y0 = 0.5*(v1[1] + v2[1]);

	r2 = (v1[0] - x0)*(v1[0] - x0) + (v1[1] - y0)*(v1[1] - y0);
}

bool Edge::Intersection( Edge* l, double* p ) {
	double denom = dx*l->dy - dy*l->dx;

	if( fabs( denom ) < 1.0e-6 ) {
		return false;
	}

	p[0] = (v1[0]*v2[1] - v1[1]*v2[0])*l->dx - dx*(l->v1[0]*l->v2[1] - l->v1[1]*l->v2[0]);
	p[1] = (v1[0]*v2[1] - v1[1]*v2[0])*l->dy - dy*(l->v1[0]*l->v2[1] - l->v1[1]*l->v2[0]);

	p[0] /= denom;
	p[1] /= denom;

	if( (p[0] - x0)*(p[0] - x0) + (p[1] - y0)*(p[1] - y0) > r2 )
		return false;

	if( (p[0] - l->x0)*(p[0] - l->x0) + (p[1] - l->y0)*(p[1] - l->y0) > r2 )
		return false;

	return true;
}

void Edge::Print() {
	cout << "[" << v1[0] << "," << v1[1] << "], [" << v2[0] << "," << v2[1] << "]" << endl;
}
