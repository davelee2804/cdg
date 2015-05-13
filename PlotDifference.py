#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import *
from pylab import *

def func( x, y ):
	theta = 0.75*np.pi
	xo = 0.5*np.cos( theta )
	yo = 0.5*np.sin( theta )
	r2 = (x-xo)*(x-xo) + (y-yo)*(y-yo)
	if np.sqrt( r2 ) < 0.4:
		return np.exp( -40.0*r2 )
	else:
		return 0.0


def eval_ij( x, y, p ):
	return p[0] + p[1]*x + p[2]*y + p[3]*x*y


nx = int(sys.argv[1])

P = np.loadtxt( 'input/Q3B2_basis_' + '%.3u'%nx + '.txt' )

XY = np.loadtxt( 'input/pgrid_' + '%.3u'%nx + '.txt' )

c = XY[:,0]
x = XY[:,1]
y = XY[:,2]

triang = Triangulation( x, y )

diff = np.zeros( len(x) )

for ii in np.arange( len(x) ):
	diff[ii] = func( x[ii], y[ii] ) - eval_ij( x[ii], y[ii], P[int(c[ii])] )

plt.tricontourf( triang, diff )
plt.colorbar()
plt.savefig( 'diff.' + '%.3u'%nx + '.png' )
plt.figure()
