class Basis {
	public:
		Basis( int _order, double* _origin );
		~Basis();
		int 		order;
		int 		nFuncs;
		double*		ci;			// coefficients
		double 		origin[2];
		double		EvalIJ( double* pt, int i );
		double 		EvalConst( double* pt );
		double 		EvalFull( double* pt );
		double		EvalWithCoeffs( double* pt, double* coeffs );
};
