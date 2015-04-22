#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#err=[0.992206,0.736548,0.536941]
#err=[0.852128,0.571278,0.335517]
err=[0.826784,0.553498,0.361314]
dx =[1.0/32,1.0/64,1.0/128]

lx = np.log(dx)
le = np.log(err)

y=np.vstack([lx,np.ones(len(lx))]).T
m=np.linalg.lstsq(y,le)[0]
print m

plt.loglog( dx, err, '-o' )
plt.show()
