class Basis {
	public:
		Basis( int _order, double* _origin, Polygon* _poly );
		~Basis();
		int 		order;
		int 		nBasisFuncs;
		double*		ci;			// coefficients
		double 		origin[2];
		Polygon*	poly;
		double		EvalIJ( double* pt, int i );
		double 		EvalConst( double* pt );
		double 		EvalFull( double* pt );
};
