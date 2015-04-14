class CDG : public CFA {
	public:
		CDG( Field* _phi, Field* _velx, Field* _vely );
		virtual ~CDG();
		double** betaInv_ij;
		virtual void Advect( double dt );
		virtual void CalcFluxes( Grid* preGrid, Field* phiTemp, double dt );
};
