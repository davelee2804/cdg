#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#Q3B2, unlimited (unnormalized)
#err_l1=[0.00333654,0.000665339,0.000172471]
#err_l2=[0.0078109,0.00146532,0.000336555]

#Q3B2 new basis
err_l1=[0.0052624,0.00094354,0.000210555]
err_l2=[0.0118542,0.00210287,0.000403864]

#Q6B3 new basis
#err_l1=[0.000332518,8.33054e-05,3.73991e-05]
#err_l2=[0.000609716,0.00016651,0.000117655]

dx =[2.0/32,2.0/64,2.0/128]

lx = np.log(dx)
le = np.log(err_l1)
y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

lx = np.log(dx)
le = np.log(err_l2)
y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

plt.loglog( dx, err_l1, '-o' )
plt.loglog( dx, err_l2, '-o' )
#plt.show()
plt.xlabel( '$\Delta x$' )
plt.ylabel( '$L_1 error$' )
plt.savefig( 'convergence_error.png' )
