class Limiter {
	public:
		Limiter( Field* _phi );
		~Limiter();
		int		order;
		Field* 	phi;
		void 	Apply();
		double 	FirstOrder( Field* phiTemp, int pi );
		double 	SecondOrder( Field* phiTemp, int pt, int dim );
};
