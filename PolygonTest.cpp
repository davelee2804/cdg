#include <iostream>
#include <cstdlib>

#include "Edge.h"
#include "Triangle.h"
#include "Basis.h"
#include "Polygon.h"

using namespace std;

#define X0 +0.0
#define Y0 +0.0

#define QUAD_ORDER 2

#define BASIS_ORDER 2

double RandX() { return ((double)6*rand())/RAND_MAX - (X0); }
double RandY() { return ((double)6*rand())/RAND_MAX - (Y0); }

int main() {
	double 		**pts;
	double		p[2]	= { -0.7, +0.3 };
	Polygon* 	poly;
	int			i;
	bool		insideTest, insideReal;
	double		nInside	= 0.0;
	int			nTests	= 1000;
	double		**pts2;
	Polygon*	poly2;
	Polygon		*poly12Int, *poly21Int;
	double		**pts3;
	Polygon*	poly3;
	Polygon		*poly13Int, *poly31Int;
	
	pts = new double*[4];
	for( i = 0; i < 4; i++ ) {
		pts[i] = new double[2];
	}
	pts[0][0] = -1.0 - (X0); pts[0][1] = -1.0 - (Y0);
	pts[1][0] = -1.0 - (X0); pts[1][1] = +1.0 - (Y0);
	pts[2][0] = +1.0 - (X0); pts[2][1] = +1.0 - (Y0);
	pts[3][0] = +1.0 - (X0); pts[3][1] = -1.0 - (Y0);

	poly = new Polygon( pts, 4, QUAD_ORDER, BASIS_ORDER );

	/* test the polygon area routine */
	cout << "area: " << poly->Area() << endl;

	/* test the point in polygon routine */
	for( i = 0; i < nTests; i++ ) {
		p[0] = RandX();
		p[1] = RandY();
		insideTest = poly->IsInside( p );
		insideReal = ( p[0] > -1.0 - (X0) && p[0] < +1 - (X0) && p[1] > -1.0 - (Y0) && p[1] < +1.0 - (Y0) ) ? true : false;
		//cout << "x: " << p->x << ", y: " << p->y << ", test: " << insideTest << ". real: " << insideReal << endl;

		if( insideTest ) {
			nInside += 1.0;
		}

		if( insideTest != insideReal ) {
			cout << "ERROR!!!" << endl;
		}
	}

	/* test the polygon intersection routine */
	pts2 = new double*[4];
	for( i = 0; i < 4; i++ ) {
		pts2[i] = new double[2];
	}
	pts2[0][0] = -1.0 - (X0) + 1.0; pts2[0][1] = -1.0 - (Y0) + 1.0;
	pts2[1][0] = -1.0 - (X0) + 1.0; pts2[1][1] = +1.0 - (Y0) + 1.0;
	pts2[2][0] = +1.0 - (X0) + 1.0; pts2[2][1] = +1.0 - (Y0) + 1.0;
	pts2[3][0] = +1.0 - (X0) + 1.0; pts2[3][1] = -1.0 - (Y0) + 1.0;
	poly2 = new Polygon( pts2, 4, QUAD_ORDER, BASIS_ORDER );
	
	poly12Int = poly->Intersection( poly2 );
	cout << "polygon 12 intersection: " << endl;
	poly12Int->Print();

	poly21Int = poly2->Intersection( poly );
	cout << "polygon 21 intersection: " << endl;
	poly21Int->Print();

	pts3 = new double*[6];
	for( i = 0; i < 6; i++ ) {
		pts3[i] = new double[2];
	}
	pts3[0][0] = +0.0; pts3[0][1] = +0.0;
	pts3[1][0] = +0.5; pts3[1][1] = +0.5;
	pts3[2][0] = +1.5; pts3[2][1] = +0.5;
	pts3[3][0] = +2.0; pts3[3][1] = +0.0;
	pts3[4][0] = +2.5; pts3[4][1] = -0.5;
	pts3[5][0] = +0.5; pts3[5][1] = -0.5;
	poly3 = new Polygon( pts3, 6, QUAD_ORDER, BASIS_ORDER );

	poly13Int = poly->Intersection( poly3 );
	cout << "polygon 13 intersection: " << endl;
	poly13Int->Print();

	poly31Int = poly3->Intersection( poly );
	cout << "polygon 31 intersection: " << endl;
	poly31Int->Print();

	delete poly2;
	delete poly3;
	delete poly12Int;
	delete poly21Int;
	delete poly13Int;
	delete poly31Int;

	/* test polygons that share a side */
	pts2[0][0] = +1.0; pts2[0][1] = +1.0;
	pts2[1][0] = +1.0; pts2[1][1] = -1.0;
	pts2[2][0] = +0.5; pts2[2][1] = -0.5;
	pts2[3][0] = +0.5; pts2[3][1] = +1.5;
	poly2 = new Polygon( pts2, 4, QUAD_ORDER, BASIS_ORDER );

	cout << "===============side share tests===============" << endl;
	cout << "polygon 2 area: " << poly2->Area() << endl;
	poly12Int = poly->Intersection( poly2 );
	cout << "polygon 12 intersection: " << endl;
	poly12Int->Print();
	cout << "polygon 12 area: " << poly12Int->Area() << endl;
	poly21Int = poly2->Intersection( poly );
	cout << "polygon 21 intersection: " << endl;
	poly21Int->Print();
	cout << "polygon 21 area: " << poly12Int->Area() << endl;

	delete poly;

	for( i = 0; i < 4; i++ ) delete[] pts[i];
	for( i = 0; i < 4; i++ ) delete[] pts2[i];
	for( i = 0; i < 6; i++ ) delete[] pts3[i];
	delete[] pts;
	delete[] pts2;
	delete[] pts3;

	return 1;
}
