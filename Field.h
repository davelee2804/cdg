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
		double		L1Error( Func* analytic, bool doNorm );
		double		L2Error( Func* analytic, bool doNorm );
		void 		Copy( Field* field );
		void		UpdateBasis();
		void 		Write( string fname, int tstep, int n );
		void 		WriteBasis( string fname, int tstep );
		void 		ReadBasis( string fname );
};
