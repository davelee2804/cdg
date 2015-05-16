class Basis {
	public:
		Basis( Polygon* _poly, int _order, double* _origin, double _dx, double _dy );
		~Basis();
		int 		order;
		int 		nFuncs;
		double		dxInv;
		double		dyInv;
		double*		ci;			// coefficients
		double*		scale;
		double*		mean;
		double 		origin[2];
		Polygon*	poly;
		void		Init();
		double		EvalIJ( double* pt, int i );
		double		EvalDerivIJ( double* pt, int i, int dim );
		double 		EvalConst( double* pt );
		double 		EvalFull( double* pt );
		double 		EvalDerivFull( double* pt, int dim );
		double		EvalWithCoeffs( double* pt, double* coeffs );
		bool		TestMean( double* volErr );
};
