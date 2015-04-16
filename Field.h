class Field {
	public:
		Field( Grid* _grid );
		~Field();
		Grid* 		grid;
		Basis**		basis;	// values are generated as an expansion of the basis
		void 		LinearInterp( double* x, double* v );
		double		IntegrateConstant(); // assumes constant cell wise values
		double		Integrate();
		void 		Copy( Field* field );
		void 		Write( string fname, int tstep );
};