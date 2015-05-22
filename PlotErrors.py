#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

q2b1_l1=[0.0103294,0.00189244,0.000376767,0.000105207]
q2b1_l2=[0.0222201,0.00436785,0.000758595,0.000199117]
dx1 =[2.0/24,2.0/48,2.0/96,2.0/192]

q4b2_l1=[0.000484247,4.96566e-05]
q4b2_l2=[0.00146153,0.000143876]
dx2 =[2.0/24,2.0/48]

lx = np.log(dx1)
le = np.log(q2b1_l1)
y=np.vstack([lx,np.ones(len(lx))]).T
m1=np.linalg.lstsq(y,le)[0]
print m1

lx = np.log(dx1)
le = np.log(q2b1_l2)
y=np.vstack([lx,np.ones(len(lx))]).T
m2=np.linalg.lstsq(y,le)[0]
print m2

lx = np.log(dx2)
le = np.log(q4b2_l1)
y=np.vstack([lx,np.ones(len(lx))]).T
m3=np.linalg.lstsq(y,le)[0]
print m3

lx = np.log(dx2)
le = np.log(q4b2_l2)
y=np.vstack([lx,np.ones(len(lx))]).T
m4=np.linalg.lstsq(y,le)[0]
print m4

plt.loglog( dx1, q2b1_l1, '-o', label='Q2B1: $L_1 error=$'+str(m1[0])[:6] )
plt.loglog( dx1, q2b1_l2, '-o', label='Q2B1: $L_2 error=$'+str(m2[0])[:6] )
plt.loglog( dx2, q4b2_l1, '-o', label='Q4B2: $L_1 error=$'+str(m3[0])[:6] )
plt.loglog( dx2, q4b2_l2, '-o', label='Q4B2: $L_2 error=$'+str(m4[0])[:6] )
plt.legend( loc='upper left' )
plt.title( 'Errors for the CDG scheme with linear basis quadratic quadrature' )
plt.xlabel( '$\Delta x$' )
plt.ylabel( 'Error' )
plt.savefig( 'convergence_error.png' )
