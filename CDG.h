typedef double ( Func ) ( double* x );

class CDG : public CFA {
	public:
		CDG( Field* _phi, Field* _velx, Field* _vely, Func* _fu, Func* _fv );
		virtual 		~CDG();
		double** 		betaInv_ij;
		void			InitBetaIJInv( Func* func );
		void 			BasisProjection( int kp, int k, double* Pij );
		void			Limiter_FirstOrder( Field* phiTemp );
		void			Limiter_SecondOrder( Field* phiTemp );
		double			PolyLimiter( Field* phiTemp, int ci );
		double			PolySecondLimiter( Field* phiTemp, int ci, int dim );
		virtual void 	Advect( double dt );
		virtual void 	CalcFluxes( Grid* preGrid, Field* phiTemp, double dt );
};
