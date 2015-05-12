#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
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

xp = np.arange(nx)*(2.0/nx) - 1.0 + 0.5*2.0/nx
yp = np.arange(nx)*(2.0/nx) - 1.0 + 0.5*2.0/nx

Xp,Yp = np.meshgrid( xp, yp )

diff = np.zeros( Xp.shape )

for ii in np.arange( Xp.shape[0] ):
	for jj in np.arange( Xp.shape[1] ):
		x = Xp[ii,jj]
		y = Yp[ii,jj]
		diff[ii,jj] = func( Xp[ii,jj], Yp[ii,jj] ) - eval_ij( Xp[ii,jj], Yp[ii,jj], P[ii*nx+jj] )

plt.contourf( Xp, Yp, diff.reshape( ( len(yp), len(xp) ) ), 100 )
plt.colorbar()
plt.savefig( 'diff' + '.' + '%.3u'%nx + '.png' )
plt.figure()
