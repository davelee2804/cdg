class Basis {
	public:
		Basis( int _order, double* _origin, double _eps );
		~Basis();
		int 	order;
		int 	nBasisFuncs;
		double*	ci;			// coefficients
		double 	origin[2];
		double	eps;
		double	EvalIJ( double* pt, int i, double* xo );
		double 	EvalConst( double* pt, double** xo );
		double 	EvalFull( double* pt, double** xo );
};
