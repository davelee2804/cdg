#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

tstep = int(sys.argv[1])
tskip = int(sys.argv[2])

xp = np.loadtxt( 'pgrid.x.txt' )
yp = np.loadtxt( 'pgrid.y.txt' )
xu = np.loadtxt( 'vgrid.x.txt' )
yu = np.loadtxt( 'vgrid.y.txt' )

velx = np.loadtxt( 'velx' + '.' + '%.4u'%0 + '.txt' )
vely = np.loadtxt( 'vely' + '.' + '%.4u'%0 + '.txt' )

velMag = np.sqrt( velx*velx + vely*vely )

Xp,Yp = np.meshgrid( xp, yp )
Xu,Yu = np.meshgrid( xu, yu )

U = velx.reshape( ( len(yu), len(xu) ) )
V = vely.reshape( ( len(yu), len(xu) ) )

plt.contourf( Xu, Yu, velMag.reshape( ( len(yu), len(xu) ) ), 100 )
plt.colorbar()
plt.quiver( Xu[::8], Yu[::8], U[::8], V[::8] )
plt.savefig( 'velMag' + '.' + '%.4u'%0 + '.png' )
plt.figure()

plt.contourf( Xu, Yu, U, 100 )
plt.colorbar()
plt.savefig( 'vel_x' + '.' + '%.4u'%0 + '.png' )
plt.figure()

plt.contourf( Xu, Yu, V, 100 )
plt.colorbar()
plt.savefig( 'vel_y' + '.' + '%.4u'%0 + '.png' )
plt.figure()

ans = np.loadtxt( 'ans' + '.' + '%.4u'%0 + '.txt' )
plt.contourf( Xp, Yp, ans.reshape( ( len(yp), len(xp) ) ), 100 )
plt.clim( 0.0, 1.0 )
plt.colorbar()
plt.savefig( 'ans' + '.' + '%.4u'%0 + '.png' )
plt.figure()

for i in np.arange( 0,tstep,tskip ):
	print 'time step: ', i
	phi = np.loadtxt( 'phi.' + '%.4u'%i + '.txt' )
	plt.contourf( Xp, Yp, phi.reshape( ( len(yp), len(xp) ) ), 100 )
	plt.clim( 0.0, 1.0 )
	plt.colorbar()
	plt.quiver( Xu[::8], Yu[::8], U[::8], V[::8] )
	plt.savefig( 'phi.' + '%.4u'%i + '.png' )
	plt.figure()

phi = np.loadtxt( 'phi.' + '%.4u'%tstep + '.txt' )
plt.contourf( Xp, Yp, phi.reshape( ( len(yp), len(xp) ) ), 100 )
plt.clim( 0.0, 1.0 )
plt.colorbar()
plt.contour( Xp, Yp, ans.reshape( ( len(yp), len(xp) ) ) )
plt.savefig( 'phi.' + '%.4u'%tstep + '.png' )
plt.figure()
