#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#Q2B2, unlimited
#err_l1=[0.0471054,0.0398877,0.084094]
#err_l2=[0.0497258,0.0507822,0.107653]

#Q2B2, unlimited (unnormalized)
#err_l1=[0.0036934,0.00312755,0.00659376]
#err_l2=[0.00985396,0.0100633,0.0213332]

#Q3B2, unlimited
#err_l1=[0.047951,0.0381468,0.0773436]
#err_l2=[0.0498745,0.0501465,0.105948]

#Q3B2, unlimited (unnormalized)
#err_l1=[0.00375985,0.00299105,0.00606446]
#err_l2=[0.00988343,0.00993733,0.0209953]

#Q3B2, unlimited (unnormalized)
err_l1=[0.00333654,0.000665339,]
#err_l2=[0.0078109,0.00146532,]

dx =[2.0/32,2.0/64,2.0/128]

lx = np.log(dx)
le = np.log(err_l1)

y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

plt.loglog( dx, err_l1, '-o' )
#plt.show()
plt.xlabel( '$\Delta x$' )
plt.ylabel( '$L_1 error$' )
plt.savefig( 'convergence_error.png' )
