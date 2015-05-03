class Grid {
	public:
		Grid( int _nx, int _ny, double _minx, double _miny, double _maxx, double _maxy, int _quadOrder, int _basisOrder, bool _internal );
		~Grid();
		int			nx;
		int			ny;
		double		minx;
		double		miny;
		double		maxx;
		double		maxy;
		int			quadOrder;
		int			basisOrder;
		bool		internal;	// true for pressure grid, false for velocity grid
		double		dx;
		double		dy;
		double		dxInv;
		double		dyInv;
		int			nVerts;
		int			nEdges;
		int			nPolys;
		double**	verts;
		Edge**		edges;
		Polygon**	polys;
		Polygon*	FindPoly( double* pt );
		int			GetPolyIndex( double* pt );
		int			EdgeCoordToIndex( int norm, int xi, int yj );
		void		EdgeIndexToCoord( int ei, int* norm, int* xi, int* yj );
		void		GetEdgePolys( int ei, Polygon* poly1, Polygon* Polygon2 );
		bool		GetEdgePolyInds( int ei, int* pinds );
		void		GetPolyEdgeInds( int pi, int* einds );
		void		GetVertPolyInds( int vi, int* cinds );
		void		GetPolyVertInds( int ci, int* vinds );
		void		UpdateEdges();
		void		UpdatePolys();
		void		UpdateTriangles();
		void		Write( string fname, int n );
};
