
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math
GeV=1.602176*10**(-10)

nop=100000000
x=np.zeros(nop)
y_gp=np.zeros(nop)
y_kl=np.zeros(nop)

l=0
test=0
for i in range (nop):
	xv=(i+1)/10.
	energy=float(GeV)
	energy=float(xv)*float(GeV)
	U = np.log(energy/ GeV * 1/200)
	x[l]=xv
	L=np.log(energy/ (1000*GeV))
	if (U >= 0 and energy >= 3 * GeV):
		cs_inel=(32.2 * (1+0.0273*U))*1e-31+32.2*0.01*pow(U,2.)*1e-31
	if (U < 0 and energy >= 3 * GeV):
		cs_inel=(32.2 * (1+0.0273*U))*1e-31
	if (energy <= 0.3 * GeV):
		cs_inel = 0
	y_gp[l]=cs_inel
	cs_inel=(34.3 + 1.88*L+0.25 *L*L)*1e-31
	y_kl[l]=cs_inel

	l=l+1

fig=plt.figure()

plt.plot(x, y_gp, color='blue')
plt.plot(x, y_kl, color='red')
plt.loglog()
plt.show()
