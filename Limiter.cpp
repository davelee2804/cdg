#include <cstdlib>
#include <cmath>

#include "Edge.h"
#include "Triangle.h"
#include "Polygon.h"
#include "Basis.h"
#include "Grid.h"
#include "Field.h"
#include "Limiter.h"

Limiter::Limiter( Field* _phi ) {
	phi = _phi;

	if( phi->basis[0]->order == 1 ) {
		order = 0;
	}
	else {
		order = ( phi->basis[0]->order == 2 ) ? 1 : 2;
	}
}

Limiter::~Limiter() {
}

void Limiter::Apply() {
	int 	pi, nx, ny, bi;
	double 	alpha_o, alpha_x, alpha_y;
	double	dx, dy, dxx, dyy, dxy;
	Basis*	basis;
	Grid*	grid	= phi->grid;

	if( !order ) {
		return;
	}

	for( pi = 0; pi < phi->grid->nPolys; pi++ ) {
		nx = pi%grid->nx;
		ny = pi/grid->nx;

		if( nx == 0 || nx == grid->nx - 1 || ny == 0 || ny == grid->ny - 1 ) {
			continue;
		}

		alpha_o = FirstOrder( pi);
		alpha_x = ( order == 2 ) ? SecondOrder( pi, 0 ) : 0.0;
		alpha_y = ( order == 2 ) ? SecondOrder( pi, 1 ) : 0.0;
		alpha_x = ( alpha_x < alpha_y ) ? alpha_x : alpha_y;
		alpha_o = ( alpha_o > alpha_x ) ? alpha_o : alpha_x;

		basis = phi->basis[pi];
		if( order == 1 && alpha_o < 1.0 - 1.0e-6 ) {
			dx = basis->ci[1];
			dy = basis->ci[basis->order];
			for( bi = 1; bi < basis->nFuncs; bi++ ) {
				basis->ci[bi] = 0.0;
			}
			basis->ci[1] = dx;
			basis->ci[basis->order] = dy;
		}
		if( order == 2 && ( alpha_o < 1.0 - 1.0e-6 || alpha_x < 1.0 - 1.0e-6 ) ) {
			dx  = basis->ci[1];
			dy  = basis->ci[basis->order];
			dxx = basis->ci[2];
			dyy = basis->ci[2*basis->order];
			dxy = basis->ci[basis->order+1];
			for( bi = 1; bi < basis->nFuncs; bi++ ) {
				basis->ci[bi] = 0.0;
			}
			basis->ci[1]              = dx;
			basis->ci[basis->order]   = dy;
			basis->ci[2]              = dxx;
			basis->ci[2*basis->order] = dyy;
			basis->ci[basis->order+1] = dxy;
		}
	}
}

double Limiter::FirstOrder( int pi ) {
	Grid*       grid            = phi->grid;
	Polygon*    poly            = grid->polys[pi];
	Basis*      basis           = phi->basis[pi];
	int         vinds[4];
	int         nx              = pi%grid->nx;
	int         ny              = pi/grid->nx;
	int         vi, ci, cj;
	double      maxPhiAtVert    = -1.0e+99;
	double      minPhiAtVert    = +1.0e+99;
	double      maxPhiInPoly    = -1.0e+99;
	double      minPhiInPoly    = +1.0e+99;
	double      phiAtVert, phiInPoly, alpha, gamma;

	grid->GetPolyVertInds( pi, vinds );

	for( vi = 0; vi < poly->n; vi++ ) {
		phiAtVert = basis->EvalFull( poly->verts[vi] ) - basis->ci[0];
		minPhiAtVert = ( phiAtVert < minPhiAtVert ) ? phiAtVert : minPhiAtVert;
		maxPhiAtVert = ( phiAtVert > maxPhiAtVert ) ? phiAtVert : maxPhiAtVert;
	}
	for( cj = ny - 1; cj < ny + 2; cj++ ) {
		for( ci = nx - 1; ci < nx + 2; ci++ ) {
			if( ci == nx && cj == ny ) {
				continue;
			}
			phiInPoly = phi->basis[cj*grid->nx+ci]->ci[0] - basis->ci[0];
			minPhiInPoly = ( phiInPoly < minPhiInPoly ) ? phiInPoly : minPhiInPoly;
			maxPhiInPoly = ( phiInPoly > maxPhiInPoly ) ? phiInPoly : maxPhiInPoly;
		}
	}

	alpha = 1.0;
	if( fabs( minPhiAtVert ) > 1.0e-6 ) {
		alpha = minPhiInPoly/minPhiAtVert > 0.0 ? minPhiInPoly/minPhiAtVert : 0.0;
	}

	gamma = 1.0;
	if( fabs( maxPhiAtVert ) > 1.0e-6 ) {
		gamma = maxPhiInPoly/maxPhiAtVert > 0.0 ? maxPhiInPoly/maxPhiAtVert : 0.0;
	}

	alpha = ( alpha < gamma ) ? alpha : gamma;
	alpha = ( alpha < 1.0 )   ? alpha : 1.0;

	return alpha;
}

double Limiter::SecondOrder( int pi, int dim ) {
	Grid*       grid            = phi->grid;
	Polygon*    poly            = grid->polys[pi];
	Basis*      basis           = phi->basis[pi];
	int         vinds[4];
	int         nx              = pi%grid->nx;
	int         ny              = pi/grid->nx;
	int         vi, ci, cj, co;
	double      maxDPhiAtVert   = -1.0e+99;
	double      minDPhiAtVert   = +1.0e+99;
	double      maxDPhiInPoly   = -1.0e+99;
	double      minDPhiInPoly   = +1.0e+99;
	double      dPhiAtVert, dPhiInPoly;
	double      alpha, gamma;

	grid->GetPolyVertInds( pi, vinds );

	for( vi = 0; vi < poly->n; vi++ ) {
		dPhiAtVert = basis->EvalDerivFull( poly->verts[vi], dim );
		minDPhiAtVert = ( dPhiAtVert < minDPhiAtVert ) ? dPhiAtVert : minDPhiAtVert;
		maxDPhiAtVert = ( dPhiAtVert > maxDPhiAtVert ) ? dPhiAtVert : maxDPhiAtVert;
	}
	for( cj = ny - 1; cj < ny + 2; cj++ ) {
		for( ci = nx - 1; ci < nx + 2; ci++ ) {
			if( ci == nx && cj == ny ) {
				continue;
			}
			co = cj*grid->nx + ci;
			dPhiInPoly = phi->basis[co]->EvalDerivFull( grid->polys[co]->origin, dim );
			minDPhiInPoly = ( dPhiInPoly < minDPhiInPoly ) ? dPhiInPoly : minDPhiInPoly;
			maxDPhiInPoly = ( dPhiInPoly > maxDPhiInPoly ) ? dPhiInPoly : maxDPhiInPoly;
		}
	}

	alpha = 1.0;
	if( fabs( minDPhiAtVert ) > 1.0e-6 ) {
		alpha = minDPhiInPoly/minDPhiAtVert > 0.0 ? minDPhiInPoly/minDPhiAtVert : 0.0;
	}

	gamma = 1.0;
	if( fabs( maxDPhiAtVert ) > 1.0e-6 ) {
		gamma = maxDPhiInPoly/maxDPhiAtVert > 0.0 ? maxDPhiInPoly/maxDPhiAtVert : 0.0;
	}

	alpha = ( alpha < gamma ) ? alpha : gamma;
	alpha = ( alpha < 1.0 )   ? alpha : 1.0;

	return alpha;
}
