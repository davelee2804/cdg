#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#gaussian ic
#q2b1_l1=[0.0103294,0.00189244,0.000376767,0.000105207]
#q2b1_l2=[0.0222201,0.00436785,0.000758595,0.000199117]
#dx1 =[2.0/24,2.0/48,2.0/96,2.0/192]
#cosine ic
q2b1_l1=[0.00360921,0.000699843,0.000147602,3.30616e-05]
q2b1_l2=[0.00900522,0.001733,0.000363742,8.21815e-05]
dx1 =[2.0/24,2.0/48,2.0/96,2.0/192]

#gaussian ic
#q4b2_l1=[0.000484247,4.96566e-05]
#q4b2_l2=[0.00146153,0.000143876]
#cosine ic
q4b2_l1=[0.000233686,3.18764e-05,4.43911e-06,9.90038e-07]
q4b2_l2=[0.00063811,9.91304e-05,1.87187e-05,4.60721e-06]
#cosine ic, half dt
q4b2_l1_0p5dt=[0.000235419,3.71071e-05,4.44266e-06,6.61061e-07]
q4b2_l2_0p5dt=[0.000641893,0.000114687,1.87089e-05,3.94163e-06]
#q4b2_l1_rk4=[0.000233715,3.66977e-05,4.44281e-06,7.08988e-07]
#q4b2_l2_rk4=[0.000638063,0.000113679,1.86816e-05,4.43494e-06]
q4b2_l1_0p5dt_rk4=[0.000235427,3.71209e-05,4.45302e-06,6.50418e-07]
q4b2_l2_0p5dt_rk4=[0.000641881,0.000114681,1.87076e-05,3.92934e-06]

dx2 =[2.0/24,2.0/48,2.0/96,2.0/192]

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

lx = np.log(dx1)
le = np.log(q4b2_l1)
y=np.vstack([lx,np.ones(len(lx))]).T
m3=np.linalg.lstsq(y,le)[0]
print m3

lx = np.log(dx1)
le = np.log(q4b2_l2)
y=np.vstack([lx,np.ones(len(lx))]).T
m4=np.linalg.lstsq(y,le)[0]
print m4

lx = np.log(dx2)
le = np.log(q4b2_l1_0p5dt_rk4)
y=np.vstack([lx,np.ones(len(lx))]).T
m5=np.linalg.lstsq(y,le)[0]
print m5

lx = np.log(dx2)
le = np.log(q4b2_l2_0p5dt_rk4)
y=np.vstack([lx,np.ones(len(lx))]).T
m6=np.linalg.lstsq(y,le)[0]
print m6

plt.loglog( dx1, q2b1_l1, '-o', label='Q2B1: $L_1 error=$'+str(m1[0])[:5] )
#plt.loglog( dx1, q2b1_l2, '-o', label='Q2B1: $L_2 error=$'+str(m2[0])[:5] )
#plt.loglog( dx1, q4b2_l1, '-o', label='Q4B2: $L_1 error=$'+str(m3[0])[:5] )
#plt.loglog( dx1, q4b2_l2, '-o', label='Q4B2: $L_2 error=$'+str(m4[0])[:5] )
plt.loglog( dx2, q4b2_l1_0p5dt_rk4, '-o', label='Q4B2: $L_1 error=$'+str(m5[0])[:5] )
#plt.loglog( dx2, q4b2_l2_0p5dt_rk4, '-o', label='Q4B2: $L_2 error (half dt)=$'+str(m6[0])[:5] )
plt.legend( loc='upper left' )
plt.title( 'Error convergence for the CDG advection scheme' )
plt.xlabel( '$\Delta x$' )
plt.ylabel( 'Error' )
plt.savefig( 'convergence_error.png' )
