class Polygon {
	public:
		Polygon( double** _verts, int _n, int quadOrder );
		~Polygon();
		int 		n;
		double 		origin[2];
		double** 	verts;
		Edge** 		edges;
		Triangle** 	tris;
		Polygon* 	Intersection( Polygon* poly );
		bool 		IsInside( double* pt );
		double 		Area();
		void 		Print();

	private:
		int			AddPoint( double** pts, int np, double* pt );
};

void MatInv( double* A, double* Ainv, int n );
void AXEB( double* A, double* x, double* b, int n );
