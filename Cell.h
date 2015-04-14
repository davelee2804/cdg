class Cell : public Polygon {
	public:
		Cell( double** verts, int _n, int _quadOrder, int _nc, double** _coords );
		virtual ~Cell();
		int			nc;
		double** 	coords;
};
