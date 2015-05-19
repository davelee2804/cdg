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

def func_dx( x, y ):
	theta = 0.75*np.pi
	xo = 0.5*np.cos( theta )
	yo = 0.5*np.sin( theta )
	r2 = (x-xo)*(x-xo) + (y-yo)*(y-yo)
	if np.sqrt( r2 ) < 0.4:
		return -40.0*2.0*(x - xo)*np.exp( -40.0*r2 )
	else:
		return 0.0

def func_dy( x, y ):
	theta = 0.75*np.pi
	xo = 0.5*np.cos( theta )
	yo = 0.5*np.sin( theta )
	r2 = (x-xo)*(x-xo) + (y-yo)*(y-yo)
	if np.sqrt( r2 ) < 0.4:
		return -40.0*2.0*(y - yo)*np.exp( -40.0*r2 )
	else:
		return 0.0

def eval_ij( x, y, p ):
	return p[0] + p[1]*x + p[2]*y + p[3]*x*y

def eval_dx( x, y, p ):
	return p[1] + p[3]*y

def eval_dy( x, y, p ):
	return p[2] + p[3]*x

nx = int(sys.argv[1])

P = np.loadtxt( 'phi_basis.0256.txt' )
dx = np.loadtxt( 'phi_deriv_0.0256.txt' )
dy = np.loadtxt( 'phi_deriv_1.0256.txt' )
XY = np.loadtxt( 'pgrid.txt' )

c = XY[:,0]
x = XY[:,1]
y = XY[:,2]

triang = Triangulation( x, y )

diff = np.zeros( len(x) )

dx_a = np.zeros( len(x) )
dy_a = np.zeros( len(x) )
dx_n = np.zeros( len(x) )
dy_n = np.zeros( len(x) )
dx_z = np.zeros( len(x) )
dy_z = np.zeros( len(x) )

for ii in np.arange( len(x) ):
	dx_n[ii] = eval_dx( x[ii], y[ii], P[int(c[ii])] )
	dy_n[ii] = eval_dy( x[ii], y[ii], P[int(c[ii])] )
	dx_a[ii] = func_dx( x[ii], y[ii] )
	dy_a[ii] = func_dy( x[ii], y[ii] )
	dx_z[ii] = dx_a[ii] - dx_n[ii]
	dy_z[ii] = dy_a[ii] - dy_n[ii]

plt.tricontourf( triang, dx, 100 )
plt.colorbar()
plt.savefig( 'dx_n.' + '%.3u'%nx + '.png' )
plt.figure()

plt.tricontourf( triang, dy, 100 )
plt.colorbar()
plt.savefig( 'dy_n.' + '%.3u'%nx + '.png' )
plt.figure()

plt.tricontourf( triang, dx_a, 100 )
plt.colorbar()
plt.savefig( 'dx_a.' + '%.3u'%nx + '.png' )
plt.figure()

plt.tricontourf( triang, dy_a, 100 )
plt.colorbar()
plt.savefig( 'dy_a.' + '%.3u'%nx + '.png' )
plt.figure()
