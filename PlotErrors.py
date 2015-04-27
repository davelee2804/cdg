#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#CFA my intersection routine
#err=[0.852128,0.571278,0.335517]

#CFA core intersection routine
#err=[0.826784,0.553498,0.361314]

#CDG constant basis funcs
err=[0.736813,0.58143,0.604775]
dx =[1.0/32,1.0/64,1.0/128]

lx = np.log(dx)
le = np.log(err)

y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

plt.loglog( dx, err, '-o' )
plt.show()
