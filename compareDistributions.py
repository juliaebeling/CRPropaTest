from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math

#~ import HadronicInteraction
h = HadronicInteraction()


nop1=100000
nop=100000-nop1/1000
x_my=np.zeros(nop)
x_e=np.zeros(nop)
x_g=np.zeros(nop)
y_my=np.zeros(nop)
y_e=np.zeros(nop)
y_g_K=np.zeros(nop)
y_g_C=np.zeros(nop)
y_c=np.zeros(nop)
y_cg=np.zeros(nop)

l=0
test=0
#~ for i in range (nop-100):
	#~ xv=(i+1)/100000.

	#~ if xv > 0.001 and xv <1:
		#~ x_my[l]=xv*TeV*10.**7./EeV
		#~ L=log(E/ TeV)
		#~ Bm= 1.75+0.204*L+0.01 * pow(L,2.)
		#~ betam=1/(1.67+0.111*L+0.0038*pow(L,2.))
		#~ km=1.07-0.086*L+0.002*pow(L,2.)
		#~ aa=(1-pow(xv,betam))/(1+km*pow(xv, betam)*(1-pow(xv,betam)))
		#~ A=Bm*log(xv)/xv*pow(aa, 4.)
		#~ B=1/log(xv)-4*betam*pow(xv,betam)/(1- pow(xv,betam))-4*km*betam*pow(xv, betam)*(1-2*pow(xv,betam))/(1+km*pow(xv,betam)*(1-pow(xv,betam)))
		#~ F=A*B
		#~ y_my[l]=F
		#~ y_c[l]=h.distribution_my1(E, xv)
		#~ l=l+1
			#~ if y_my[l]==y_c[l]:
				#~ test=test+1

#~ for j in range (5):
	#~ l=0
	#~ E=10**(j-1)*TeV
	#~ for i in range (nop-1000):
		#~ xv=(i+1)/100000.
		#~ if xv > 0.001 and xv <1:
			#~ x_g[l]=xv*TeV
			#~ y_g_K[l]=h.distribution_gamma(E, xv)
			#~ y_cg[l]=2*h.distribution_Carceller_g(E, xv, 1., 5.8, 5.3)
			#~ l=l+1
	#~ integral_K=np.sum(y_g_K)/len(x_g)*(np.max(x_g)-np.min(x_g))
	#~ y_g_K=y_g_K/integral_K
	#~ integral_C=np.sum(y_cg)/len(x_g)*(np.max(x_g)-np.min(x_g))
	#~ y_cg=y_cg/integral_C

l=0
for i in range (nop):
	xv=(i+1)/100000.
	x_g[l]=xv
	y_g_K[l]=h.CrossSection_Kelner(xv*PeV)
	y_cg[l]=h.CrossSection_Carceller(xv*PeV)
	y_c[l]=h.CrossSection_Galprop(xv*PeV)
	l=l+1


fig=plt.figure()
plt.xlabel('Energy E/TeV')
plt.ylabel('Cross Section in mb')
plt.plot(x_g, y_cg, color='blue', label='Carceller')
plt.plot(x_g, y_g_K, color='red', label='Kelner')
plt.plot(x_g, y_c, color='green', label='Galprop')
plt.legend(loc='upper right')
plt.xscale('log')
plt.title ('Inelastic Cross Sections based on Different Models')
fig.savefig('/home/home1/jce/KelnerCarceller/CS.png')
plt.close(fig)
#~ plt.show()
