class Basis {
	public:
		Basis( Polygon* _poly, int _order, double* _origin, double _dx, double _dy );
		~Basis();
		int 		order;
		int 		nFuncs;
		double		dxInv;
		double		dyInv;
		double		aInv;
		double*		ci;			// coefficients
		double 		origin[2];
		Polygon*	poly;
		double		EvalIJ( double* pt, int i );
		double 		EvalConst( double* pt );
		double 		EvalFull( double* pt );
		double		EvalWithCoeffs( double* pt, double* coeffs );
};
