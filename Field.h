typedef double ( Func ) ( double* x );

class Field {
	public:
		Field( Grid* _grid );
		~Field();
		Grid* 		grid;
		Basis**		basis;	// values are generated as an expansion of the basis
		double		EvalAtCoord( double* x );
		void 		LinearInterp( double* x, double* v );
		double		IntegrateConstant(); // assumes constant cell wise values
		double		Integrate();
		double		L1Error( Func* analytic );
		double		L2Error( Func* analytic );
		void 		Copy( Field* field );
		void		UpdateBasis();
		void 		Write( string fname, int tstep, int n );
};
