#include <iostream>

#include "Edge.h"

using namespace std;

int main() {
	Edge* 	l1	= new Edge( -1.0, -1.0, +1.0, +1.0 );
	Edge*	l2	= new Edge( -1.0, +1.0, +1.0, -1.0 );
	Edge*	l3	= new Edge( +0.01, +0.01, +2.0, +2.0 );
	Point* 	p	= new Point();

	l1->Print();
	l2->Print();
	l3->Print();

	//cout << "[l1, l2] intersection: " << l1->Intersection( l2, p ) << endl;
	if( l1->Intersection( l2, p ) )
		cout << p->x << "," << p->y << endl;
	cout << "[l2, l3] intersection: " << l2->Intersection( l3, p ) << endl;

	delete l1;
	delete l2;
	delete l3;
	//delete p;
}
