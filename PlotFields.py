#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import *
from pylab import *

def analytic( x, y ):
	xo = 0.5*np.cos( 0.75*np.pi )
	yo = 0.5*np.sin( 0.75*np.pi )
	r2 = (x - xo)*(x - xo) + (y - yo)*(y - yo)
	if( np.sqrt( r2 ) ) < 0.4:
		return 0.5*( 1.0 + np.cos( np.pi*np.sqrt(r2)/0.4 ) )
	else:
		return 0.0

def basis_eval_linear( x, y, p, nx ):
	#dx = 4.0/nx
	#return p[0] + p[1]*x/dx + p[2]*y/dx
	return p[0] + p[1]*x + p[2]*y

def basis_eval_quadratic( x, y, p, nx ):
	#dx = 4.0/nx
	#return p[0] + p[1]*x/dx + p[2]*y/dx + p[3]*x*x/dx/dx/2 + p[4]*x*y/dx/dx/2 + p[5]*y*y/dx/dx
	return p[0] + p[1]*x + p[2]*y + p[3]*x*x + p[4]*x*y + p[5]*y*y

tstep = int(sys.argv[1])
tskip = int(sys.argv[2])
nx    = int(sys.argv[3])

XY = np.loadtxt( 'pgrid.txt' )
p = XY[:,0]
x = XY[:,1]
y = XY[:,2]

triang = Triangulation( x, y )

z = np.zeros( len(x) )
a = np.zeros( len(x) )

for i in np.arange( 0, tstep, tskip ):
	print 'time step: ', i
	phi = np.loadtxt( 'phi_basis.' + '%.4u'%i + '.txt' )
	for j in np.arange( len(x) ):
		z[j] = basis_eval_quadratic( x[j], y[j], phi[int(p[j])], nx )

	plt.tricontourf( x, y, z, 100 )
	plt.clim( 0.0, 1.0 )
	plt.colorbar()
	plt.savefig( 'phi.' + '%.4u'%i + '.png' )
	plt.figure()

for i in np.arange( len(x) ):
	a[i] = analytic( x[i], y[i] )

phi = np.loadtxt( 'phi_basis.' + '%.4u'%tstep + '.txt' )
for j in np.arange( len(x) ):
	z[j] = basis_eval_quadratic( x[j], y[j], phi[int(p[j])], nx )

plt.tricontourf( x, y, z, 100 )
plt.clim( 0.0, 1.0 )
plt.colorbar()
plt.tricontour( x, y, a )
plt.savefig( 'phi.' + '%.4u'%tstep + '.png' )
plt.figure()


