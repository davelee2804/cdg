#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#Q3B2 unlimited (unnormalized) (with cross terms)
#err_l1=[0.00333654,0.000665339,0.000172471]
#err_l2=[0.0078109,0.00146532,0.000336555]

#Q3B2 unlimited (unnormalized)
#err_l1=[0.0052624,0.00094354,0.000210555]
#err_l2=[0.0118542,0.00210287,0.000403864]
#dx =[2.0/32,2.0/64,2.0/128]

err_l1=[0.0104051,0.00190689,0.000376642,0.000103886]
err_l2=[0.0222779,0.00439083,0.000767012,0.000201618]
dx =[2.0/24,2.0/48,2.0/96,2.0/192]

#Q6B3 new basis
#err_l1=[0.000332518,8.33054e-05,3.73991e-05]
#err_l2=[0.000609716,0.00016651,0.000117655]

lx = np.log(dx)
le = np.log(err_l1)
y=np.vstack([lx,np.ones(len(lx))]).T
m1=np.linalg.lstsq(y,le)[0]
print m1

lx = np.log(dx)
le = np.log(err_l2)
y=np.vstack([lx,np.ones(len(lx))]).T
m2=np.linalg.lstsq(y,le)[0]
print m2

plt.loglog( dx, err_l1, '-o', label='$L_1 error=$'+str(m1[0]) )
plt.loglog( dx, err_l2, '-o', label='$L_2 error=$'+str(m2[0]) )
plt.legend( loc='upper left' )
plt.title( 'Errors for the CDG scheme with linear basis quadratic quadrature' )
plt.xlabel( '$\Delta x$' )
plt.ylabel( 'Error' )
plt.savefig( 'convergence_error.png' )
