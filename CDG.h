class CDG : public CFA {
	public:
		CDG( Field* _phi, Field* _velx, Field* _vely );
		virtual 		~CDG();
		double** 		betaInv_ij;
		double			phiMax;
		double			phiMin;
		void 			BasisProjection( int kp, int k, double* Pij );
		void			Limiter( Field* phiTemp );
		virtual void 	Advect( double dt );
		virtual void 	CalcFluxes( Grid* preGrid, Field* phiTemp, double dt );
};
