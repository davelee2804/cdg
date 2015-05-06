#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#CFA my intersection routine
#err=[0.852128,0.571278,0.335517]

#CFA core intersection routine
#err=[0.826784,0.553498,0.361314]

#CDG second order (bill's limiter)
#err=[0.679821,0.498348,0.914343]

#CDG second order (no limiter)
#err=[0.286685,0.179829,0.178711]

#CDG second order (kuzmin limiter)
err=[0.680037,0.498419,0.914583]

#CDG second order (venkatakrishnan limiter)
#err=[1.38352,0.706516,0.962705]

dx =[1.0/32,1.0/64,1.0/128]

lx = np.log(dx)
le = np.log(err)

y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

plt.loglog( dx, err, '-o' )
#plt.show()
plt.xlabel( '$\Delta x$' )
plt.ylabel( '$L_1 error$' )
plt.savefig( 'convergence_error.png' )
