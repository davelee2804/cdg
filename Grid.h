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
		int			nCells;
		double**	verts;
		Edge**		edges;
		Cell**		cells;
		Cell*		FindCell( double* pt );
		int			GetCellIndex( double* pt );
		int			EdgeCoordToIndex( int norm, int xi, int yj );
		void		EdgeIndexToCoord( int ei, int* norm, int* xi, int* yj );
		void		GetEdgeCells( int ei, Cell* cell1, Cell* Cell2 );
		bool		GetEdgeCellInds( int ei, int* pinds );
		void		GetCellEdgeInds( int pi, int* einds );
		void		GetVertCellInds( int vi, int* cinds );
		void		GetCellVertInds( int ci, int* vinds );
		void		UpdateEdges();
		void		UpdateCells();
		void		Write( string fname );
};
