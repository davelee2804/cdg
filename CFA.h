typedef double ( Func ) ( double* x );

class CFA {
	public:
		CFA( Field* _phi, Field* _velx, Field* _vely, Func* _fu, Func* _fv );
		~CFA();
		Field* 		phi;
		Field* 		velx;
		Field* 		vely;
		Func*		fu;
		Func*		fv;
		double** 	pts;
		void 		Advect( double dt );
		void 		CalcChars( Grid* preGrid, double dt );
		void 		CalcFluxes( Grid* preGrid, Field* phiTemp, double dt );
		void 		TraceEuler( double dt, int dir, double* xi, double* xf );
		void 		TraceRK2( double dt, int dir, double* xi, double* xf );
		double 		GetNorm( double* a, double* b, double* c );
		void		CheckBounds( double* pt );
		Polygon* 	CreatePreImage( int ei, Grid* grid, Grid* preGrid, int* into, int* from, int* pinds );
		Polygon* 	Intersection( Polygon* poly1, Polygon* poly2 );
};
