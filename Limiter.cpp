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
	int 	pi, bi;
	double 	alpha_o, alpha_x, alpha_y;
	double	dx, dy, dxx, dyy, dxy;
	Basis*	basis;

	if( !order ) {
		return;
	}

	for( pi = 0; pi < phi->grid->nPolys; pi++ ) {
		alpha_o = FirstOrder( pi );
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
	Grid* 		grid 	= phi->grid;
	Polygon*	poly	= grid->polys[pi];
	Basis*		basis	= phi->basis[pi];
	int 		vinds[4];
	int			nx		= pi%grid->nx;
	int			ny		= pi/grid->nx;
	int			xmin	= ( nx > 0 ) ? nx - 1  : nx;
	int			ymin	= ( ny > 0 ) ? ny - 1  : ny;
	int			xmax	= ( nx < grid->nx - 1 ) ? nx + 1 : nx;
	int			ymax	= ( ny < grid->ny - 1 ) ? ny + 1 : ny;
	int			ci, cj, co, vi;
	double		phiMin	= +1.0e+99;
	double		phiMax	= -1.0e+99;
	double		alpha	= 1.0;
	double		alpha_i;
	double		phiOrig, phiVert, ratio;

	for( cj = ymin; cj <= ymax; cj++ ) {
		for( ci = xmin; ci <= xmax; ci++ ) {
			co = cj*grid->nx + ci;
			phiOrig = phi->basis[co]->ci[0];
			phiMax = ( phiOrig > phiMax ) ? phiOrig : phiMax;
			phiMin = ( phiOrig < phiMin ) ? phiOrig : phiMin;
		}
	}

	grid->GetPolyVertInds( pi, vinds );

	phiOrig = basis->ci[0];
	for( vi = 0; vi < 4; vi++ ) {
		phiVert = basis->EvalFull( poly->verts[vi] );
		alpha_i = 1.0;
		if( phiVert > phiOrig + 1.0e-6 ) {
			ratio = (phiMax - phiOrig)/(phiVert - phiOrig);
			alpha_i = ( 1.0 < ratio ) ? 1.0 : ratio;
		}
		else if( phiVert < phiOrig - 1.0e-6 ) {
			ratio = (phiMin - phiOrig)/(phiVert - phiOrig);
			alpha_i = ( 1.0 < ratio ) ? 1.0 : ratio;
		}
		alpha = ( alpha_i < alpha ) ? alpha_i : alpha;
	}

	return alpha;
}

double Limiter::SecondOrder( int pi, int dim ) {
	Grid* 		grid 	= phi->grid;
	Polygon*	poly	= grid->polys[pi];
	Basis*		basis	= phi->basis[pi];
	int 		vinds[4];
	int			nx		= pi%grid->nx;
	int			ny		= pi/grid->nx;
	int			xmin	= ( nx > 0 ) ? nx - 1  : nx;
	int			ymin	= ( ny > 0 ) ? ny - 1  : ny;
	int			xmax	= ( nx < grid->nx - 1 ) ? nx + 1 : nx;
	int			ymax	= ( ny < grid->ny - 1 ) ? ny + 1 : ny;
	int			ci, cj, co, vi;
	double		dPhiMin	= +1.0e+99;
	double		dPhiMax	= -1.0e+99;
	double		alpha	= 1.0;
	double		alpha_i;
	double		dPhi, dPhiVert, ratio;

	for( cj = ymin; cj <= ymax; cj++ ) {
		for( ci = xmin; ci <= xmax; ci++ ) {
			co = cj*grid->nx + ci;
			dPhi = phi->basis[co]->EvalDerivFull( grid->polys[co]->origin, dim );
			dPhiMax = ( dPhi > dPhiMax ) ? dPhi : dPhiMax;
			dPhiMin = ( dPhi < dPhiMin ) ? dPhi : dPhiMin;
		}
	}

	grid->GetPolyVertInds( pi, vinds );

	dPhi = basis->EvalDerivFull( poly->origin, dim );
	for( vi = 0; vi < 4; vi++ ) {
		dPhiVert = basis->EvalDerivFull( poly->verts[vi], dim );
		alpha_i = 1.0;
		if( dPhiVert > dPhi + 1.0e-6 ) {
			ratio = (dPhiMax - dPhi)/(dPhiVert - dPhi);
			alpha_i = ( 1.0 < ratio ) ? 1.0 : ratio;
		}
		else if( dPhiVert < dPhi - 1.0e-6 ) {
			ratio = (dPhiMin - dPhi)/(dPhiVert - dPhi);
			alpha_i = ( 1.0 < ratio ) ? 1.0 : ratio;
		}
		alpha = ( alpha_i < alpha ) ? alpha_i : alpha;
	}

	return alpha;
}
