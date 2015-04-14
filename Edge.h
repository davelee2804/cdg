#include <iostream>

using namespace std;

class Edge {
	public:
		Edge( double* x1, double* x2 );
		~Edge();
		double GetY( double x );
		double GetX( double y );
		bool Intersection( Edge* l, double* p );
		void Print();
		double m;
		double c;
		double dx;
		double dy;
		double x0;
		double y0;
		double r2;
		double v1[2];
		double v2[2];
};
