class Limiter {
	public:
		Limiter( Field* _phi );
		~Limiter();
		int		order;
		Field* 	phi;
		void 	Apply();
		double 	FirstOrder( int pi );
		double 	SecondOrder( int pt, int dim );
};
