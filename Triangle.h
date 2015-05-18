class Triangle {
	public:
		Triangle( double* x1, double* x2, double* x3, int _order );
		~Triangle();
		double 		a[2];
		double 		b[2];
		double 		c[2];
		double		area;
		Edge* 		ab;
		Edge* 		bc;
		Edge* 		ca;
		int 		order;
		int 		nQuadPts;
		double* 	wi;		// quadrature weights
		double** 	qi;		// quadrature points
		void		Init();
		void 		Quadrature();
		double 		CrossProduct( int origin );
};
