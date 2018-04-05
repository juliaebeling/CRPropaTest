from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math

#~ import HadronicInteraction
h = HadronicInteraction()


nop1=10000000000
nop=100000
x_my=np.zeros(nop)
x_e=np.zeros(nop)
x_g=np.zeros(nop)
y_my=np.zeros(nop)
y_e=np.zeros(nop)
y_g_K=np.zeros(nop)
y_g_K1=np.zeros(nop)
y_g_K2=np.zeros(nop)
y_g_K3=np.zeros(nop)
y_g_C=np.zeros(nop)
y_g_C1=np.zeros(nop)
y_g_C2=np.zeros(nop)
y_c=np.zeros(nop)
y_cg=np.zeros(nop)
E=10*TeV
l=0
test=0
for i in range (nop-100):
    xv=(i+1)/100000.
    if (xv > (0.001) and xv <1):
        x_my[l]=xv
        L=log(E/ TeV)
        Bm= 1.75+0.204*L+0.01 * pow(L,2.)
        betam=1/(1.67+0.111*L+0.0038*pow(L,2.))
        km=1.07-0.086*L+0.002*pow(L,2.)
        y=xv/0.427
        aa=(1-pow(y,betam))/(1+km*pow(y, betam)*(1-pow(y,betam)))
        A=Bm*log(y)/y*pow(aa, 4.)
        B=1/log(y)-4*betam*pow(y,betam)/(1- pow(y,betam))-4*km*betam*pow(y, betam)*(1-2*pow(y,betam))/(1+km*pow(y,betam)*(1-pow(y,betam)))
        F=A*B
        
        if xv < 0.427:
            y_c[l]=h.distribution_my1(E, xv)
            y_my[l]=h.distribution_e(E,xv)
            y_g_C[l]=h.distribution_my1(E, xv)
        y_g_K[l]=(h.distribution_gamma(E, xv))
        y_g_C[l]=h.distribution_e(E,xv)+y_g_C[l]
        l=l+1
integral_K=np.sum(y_g_K)/len(x_my)*(np.max(x_my)-np.min(x_my))
        
print integral_K
integral_my=np.sum(y_my)/len(x_my)*(np.max(x_my)-np.min(x_my))
print integral_my
print h.number_my1(0.1*TeV)



#~ for i in range (1):
	#~ l=0
	#~ E=TeV*10**3
	#~ for i in range (nop):
		#~ xv=(i+1)*10**(-9)
		#if xv > 0.0001 and xv <1:
		#x_g[l]=xv
			#~ y_g_K[l]=h.distribution_gamma(E, xv)*xv**2
			#~ if xv < 0.427:
				#~ y_g_K1[l]=h.distribution_my1(E, xv)*xv**2
			#~ y_g_K3[l]=(h.distribution_e(E, xv)*xv**2+y_g_K1[l])
			#~ y_g_K2[l]=h.distribution_e(E, xv)*xv**2

	#~ integral_K=np.sum(y_g_K)/len(x_g)*(np.max(x_g)-np.min(x_g))
	#~ y_g_K=y_g_K/integral_K
	#~ integral_K1=np.sum(y_g_K1)/len(x_g)*(np.max(x_g)-np.min(x_g))
	#~ y_g_K1=y_g_K1/integral_K1
	#~ integral_K2=np.sum(y_g_K2)/len(x_g)*(np.max(x_g)-np.min(x_g))
	#~ y_g_K2=y_g_K2/integral_K2

		#~ y_g_C[l]=h.distribution_Carceller(E, xv, 1, 15.5, 8.0)*xv
		#~ y_g_C1[l]=h.distribution_Carceller(E, xv, 4, 14.8, 7.3)*xv
		#~ y_g_C2[l]=h.distribution_Carceller(E, xv, 56, 20.0, 5.8)*xv
		#~ l=l+1
	#l=0
#~ #Cross Section models
#for i in range (nop):
#    xv=(i+1)/100000.
#    x_g[l]=xv*1000
#    y_g_K[l]=h.CrossSection_Kelner(xv*PeV)*10**31
#    y_cg[l]=h.CrossSection_Carceller(xv*PeV)*10**31
#    y_c[l]=h.CrossSection_Galprop(xv*PeV)*10**31
#    l=l+1



	#~ #Cross Section pion charge
	#~ for i in range (nop):
	#~ xv=(i+1)/10000000.
	#~ x_g[i]=xv*1000000000
	#~ xv=xv*1000000000*TeV
	#~ y_g_K[i]=1/(0.007+0.1*np.log(xv)/xv+0.3/xv**2)
	#~ y_cg[i]=1/(0.00717+0.0652*np.log(xv)/xv+0.162/xv**2)
	#~ y_c[i]=1/(0.00456+0.0846/xv**0.5+0.577/xv**1.5)




fig=plt.figure()
plt.xlim(10**(-3), 1)
plt.ylim(0.001, 0.1)
#plt.xlabel('Energy E/E_p')
#plt.ylabel('f*x**2')
#plt.plot(x_g, y_g_K, color='blue', label='$\gamma$ Kelner et al 2006')
plt.plot(x_my, y_g_C*x_my**2, color='dodgerblue', label='Hydrogen')
#plt.plot(x_g, y_g_K2, color='green', label='   Kelner et al 2006')
#plt.plot(x_g, y_g_C2, color='green', label='Helium')
#plt.plot(x_g, y_g_K3, color='red', label='    Kelner et al 2006')
#plt.plot(x_g, y_g_C1, color='red', label='Iron')

plt.plot(x_my, y_g_K*x_my**2, label='Gamma', color='blue')
plt.plot(x_my, y_my*x_my**2, label='Myp', color='red')
plt.plot(x_my, y_c*x_my**2, label='MyHI', color='green')
#plt.plot(Sb_x, Sb_y, marker='.', linewidth=0, color='grey', label='SIBYLL')
#plt.plot(R_x, R_y, marker='.', linewidth=0, color='black', label='IHEP')


plt.legend(loc='upper left')
plt.loglog()
plt.title ('Inelastic Cross Section According to Different Models')
#fig.savefig('/home/home1/jce/KelnerCarceller/CSdata.png')
plt.show()
#plt.close(fig)

Sb_x=np.array([0.05064603544,
               0.08111308308,
               0.1299081397,
               0.2080567538,
               0.3332170941,
               0.5198871593,
               0.8326347754,
               1.299081397,
               2.080567538,
               3.420510352,
               5.198871593,
               8.326347754,
               13.33521432,
               20.80567538,
               33.32170941,
               53.36699231,
               83.26347754,
               133.3521432,
               213.5725606,
               333.2170941,
               533.6699231,
               854.7088126])

Sb_y=np.array([30.67567568,
               31.08108108,
               31.48648649,
               32.02702703,
               32.56756757,
               33.37837838,
               34.18918919,
               35,
               35.94594595,
               37.02702703,
               38.10810811,
               39.18918919,
               40.40540541,
               41.75675676,
               43.10810811,
               44.45945946,
               46.08108108,
               47.83783784,
               50.54054054,
               52.02702703,
               54.59459459,
               57.02702703])

R_x=np.array([0.00877368056,
              0.01232846739,
              0.01201006767,
              0.023101297,
              0.03080607466,
              0.05772495943,
              0.1,
              0.1687612476,
              0.2026834004,
              0.2848035868,
              0.3898603703,
              0.4933803153,
              1.026511068,
              1.480657277,
              2.026834004])

R_y=np.array([30,
              29.72972973,
              28.91891892,
              30.54054054,
              30,
              31.62162162,
              31.75675676,
              31.62162162,
              32.2972973,
              32.16216216,
              32.7027027,
              32.97297297,
              34.59459459,
              35.13513514,
              35.54054054])
