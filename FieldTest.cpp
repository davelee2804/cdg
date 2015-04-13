#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Basis.h"
#include "Polygon.h"
#include "Grid.h"
#include "Field.h"

using namespace std;

#define QUAD_ORDER 2
#define BASIS_ORDER 2

int main() {
	int 	nx		= 2;
	int		ny		= 2;
	double	v, x[2];
	int 	i;
	Grid*	grid 	= new Grid( nx, ny, -1.0, -1.0, +1.0, +1.0, QUAD_ORDER, BASIS_ORDER );
	Field*	field	= new Field( grid, 1, false );

	//field->vals[(ny/2)*(nx+1)+(nx/2)][0] = 1.0;
	for( i = 0; i < grid->nPts; i++ ) {
		field->vals[i][0] = grid->verts[i][0]*grid->verts[i][1];
	}

	x[0] = -0.5;
	x[1] = +0.5;

	field->Interp( x, &v );
	cout << "test val: " << v << endl;

	delete grid;
	delete field;

	return 1;
}
