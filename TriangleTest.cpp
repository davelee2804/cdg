#include <iostream>

#include "Edge.h"
#include "Triangle.h"

using namespace std;

#define TRI_ORDER 2

int main() {
	int 		i;
	double		a[2]	= { 0.0, 0.0 };
	double		b[2]	= { 0.0, 2.0 };
	double		c[2]	= { 2.0, 0.0 };
	Triangle* 	tri 	= new Triangle( a, b, c, TRI_ORDER );

	cout << "triangle area: " << tri->Area() << endl;

	cout << "triangle vertices: " << endl;
	cout << "[" << tri->a[0] << "," << tri->a[1] << "]" << endl;
	cout << "[" << tri->b[0] << "," << tri->b[1] << "]" << endl;
	cout << "[" << tri->c[0] << "," << tri->c[1] << "]" << endl;

	cout << "quadrature coords: " << endl;
	for( i = 0; i < tri->nq; i++ ) {
		cout << "[" << tri->qi[i][0] << "," << tri->qi[i][1] << "]" << endl;
	}

	return 1;
}
