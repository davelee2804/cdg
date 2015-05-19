#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#Q3B2, unlimited (unnormalized)
#err_l1=[0.00333654,0.000665339,0.000172471]
err_l2=[0.0078109,0.00146532,0.000336555]

dx =[2.0/32,2.0/64,2.0/128]

lx = np.log(dx)
le = np.log(err_l2)

y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

plt.loglog( dx, err_l2, '-o' )
#plt.show()
plt.xlabel( '$\Delta x$' )
plt.ylabel( '$L_1 error$' )
plt.savefig( 'convergence_error.png' )
