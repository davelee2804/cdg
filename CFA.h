class CFA {
	public:
		CFA( Field* _phi, Field* _velx, Field* _vely );
		~CFA();
		Field* phi;
		Field* velx;
		Field* vely;
		void Advect( double dt );
		void CalcChars( Grid* preGrid, double dt );
		void CalcFluxes( Grid* preGrid, Field* phiTemp, double dt );
		void TraceRK2( double dt, int dir, double* xi, double* xf );
		double GetNorm( double* a, double* b, double* c );
};
