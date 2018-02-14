from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math



#Number of particles
particles = 100000000

#Steplength
steplength = 10000000000000

#Candidate
c = Candidate()


#Energy
E = 1 * TeV
Energyvector=[]
source= Source()
source.add(SourceEnergy(E))
source.add(SourceParticleType(nucleusId(1,1)))
source.add(SourcePosition(Vector3d(0)))
output = TextOutput("test2.txt")
output.disableAll()
output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
obs=Observer()

obs.add(ObserverTimeEvolution(0, 1000 * kpc, 100000))
obs.onDetection(output)

m = ModuleList()

m.add(SimplePropagation(steplength * kpc,steplength* kpc))
m.add(HadronicInteraction())
m.add(obs)
m.setShowProgress(True)
m.run(source, particles, True)

e_dist_my=[]
e_dist_e=[]
e_dist_ne=[]
e_dist_g=[]
e_dist_p=np.zeros(1000)
e_dist_x=np.zeros(1000)
data=[]
data = pd.read_csv('test2.txt',
       names=['SN','ID','E','SNO','SN1'], delimiter='\t', comment='#',
       usecols=["SN", "ID", "E"])
#~ data.sort_values('E', ascending=True)
data['ID']=data.ID
data['E']=data.E

data.ID.tolist()
data.E.tolist()




l=0
x_test=[]
#Einlesen der verschiedenen Teilchen
for i in range (particles):
	if data.ID[i]==11:
		e_dist_e.append(data.E[i]*EeV)
	if data.ID[i]==14:
		e_dist_my.append(data.E[i]*EeV)
	if data.ID[i]==12:
		e_dist_ne.append(data.E[i]*EeV)
	if data.ID[i]==22:
		e_dist_g.append(data.E[i]*EeV)
		l=l+1
		
edist_min_g=np.min(e_dist_g)
edist_max_e=math.ceil(np.max(e_dist_e))
edist_max_my=math.ceil(np.max(e_dist_my))
edist_max_ne=math.ceil(np.max(e_dist_ne))/1.00000
#~ for i in range(len(e_dist_p)):
	#~ print e_dist_ne[i]

nop1=100000
nop=100000-nop1/1000
x=np.zeros(nop)
x_my=np.zeros(nop)
x_g=np.zeros(nop)
y=np.zeros(nop)
y_my=np.zeros(nop)
y_g=np.zeros(nop)
x_compare=[]
x_compares=[]
y_compare=[]
y_sim=[]
#Verteilungsfunktion Elektron
l=0
n=1
for i in range (nop):
	xv=(i+1)/100000.
	if xv >= 0.001 and xv <1:
		x[l]=xv*E*1.
		L=log(E/ TeV)
		Be= 1/(69.5+2.65*L+0.3*L**2)
		betae=1/(0.201+0.062*L+0.00042*L**2)**0.25
		ke=(0.279 + 0.141 *L + 0.0172* L**2)/(0.3+ (2.3+L)**2.)

		F=Be*(1+ke*log(xv)**2.)**3. /(xv*(1+0.3/xv**betae))*(-log(xv)**5.)
		y[l]=F
		l=l+1
	#~ F=xv**4.
	
nob=100
e_test=np.max(e_dist_e)

z= np.histogram(e_dist_e, bins =np.arange(0, edist_max_e, edist_max_e/nob), normed=True, weights=None, density=None)

integral=np.sum(y)/(len(x)*0.1)
y=y/integral


#~ for i in range (1, nob-1):
	#~ nob=nob*1.
	#~ xs=i/nob*1.
	#~ xa=i/nob*1.*(np.max(x)/np.max(edist_max_e))**(-1)
	#~ L=log(E/ TeV)
	#~ Be= 1/(69.5+2.65*L+0.3*L**2)
	#~ betae=1/(0.201+0.062*L+0.00042*L**2)**0.25
	#~ ke=(0.279 + 0.141 *L + 0.0172* L**2)/(0.3+ (2.3+L)**2.)
	#~ F=Be*(1+ke*log(xa)**2.)**3. /(xa*(1+0.3/xa**betae))*(-log(xa)**5.)
	#~ yc=F/integral*1.
	#~ x_compares.append(xs*edist_max_e)
	#~ x_compare.append(xa*TeV*10.**7./EeV)
	#~ y_compare.append(yc)
	#~ y_sim.append(z[0][i-1])
	
l=0
for i in range (nop):
	xv=(i+1)/100000.
	if xv >= 0.001 and xv <1:
		x_my[l]=xv*E
		L=log(E/ TeV)
		Bm= 1.75+0.204*L+0.01 * pow(L,2.)
		betam=1/(1.67+0.111*L+0.0038*pow(L,2.))
		km=1.07-0.086*L+0.002*pow(L,2.)
		aa=(1-pow(xv,betam))/(1+km*pow(xv, betam)*(1-pow(xv,betam)))
		A=Bm*log(xv)/xv*pow(aa, 4.)
		B=1/log(xv)-4*betam*pow(xv,betam)/(1- pow(xv,betam))-4*km*betam*pow(xv, betam)*(1-2*pow(xv,betam))/(1+km*pow(xv,betam)*(1-pow(xv,betam)))
		F=A*B
		y_my[l]=F
		l=l+1

integral=np.sum(y_my)/len(x_my)*(np.max(x_my)-np.min(x_my))
y_my=y_my/integral

l=0
for i in range (nop):
	xv=(i+1)/100000.
	if xv >= 0.001 and xv <1:
		x_g[l]=xv*E
		L=log(E / TeV)
		Bg=1.3+0.14*L+0.011*L**2
		betag=1/(1.79+0.11*L+0.008*L**2)
		kg=1/(0.801+0.049*L+0.014*L**2)
		A=Bg*np.log(xv)/xv
		B=(1-xv**betag)/(1+kg*xv**betag*(1-xv*betag))
		C=1/log(xv)-4*betag*xv**betag/(1-xv**betag)-4*kg*betag*xv**betag*(1-2*xv**betag)/(1+kg*xv**betag*(1-xv**betag))
		F=A*B**4*C
		y_g[l]=F
		l=l+1

integral=np.sum(y_g)/len(x_g)*(np.max(x_g)-np.min(x_g))
y_g=y_g/integral

fig=plt.figure()

#~ #plt.hist(e_dist_my, bins = np.arange(0, edist_max_my, edist_max_my/20 ), color='blue')
#~ plt.hist(e_dist_e, bins = np.arange(0, edist_max_e, edist_max_e/100 ), normed=True, color='green')
Edges = np.arange(0, edist_max_e, edist_max_e/100 )
Edges = np.logspace(-10, -6, 1000)
h, binEdges_e = np.histogram(e_dist_e, bins=Edges, normed=True)
my, binEdges_my = np.histogram(e_dist_my, bins=Edges, normed=True)
g, binEdges_g=np.histogram(e_dist_g, bins=Edges, normed=True)
binMiddle_e = (binEdges_e[1:]-binEdges_e[:-1])/2.+binEdges_e[:-1]
binMiddle_my = (binEdges_my[1:]-binEdges_my[:-1])/2.+binEdges_my[:-1]
binMiddle_g=(binEdges_g[1:]-binEdges_g[:-1])/2.+binEdges_g[:-1]
#~ plt.plot(binMiddle_e, h, linewidth=0., marker='s')
#~ plt.plot(binMiddle_my, my, linewidth=0., marker='s')
#~ print g
plt.plot(binMiddle_g, g, linewidth=0., marker='s')
#~ plt.plot(binRight, h, linewidth=0., marker='s')
#~ plt.scatter(x_compares, y_sim)
#~ plt.plot(x_my, y_my, color='red')
plt.plot(x_g, y_g, color='red')
plt.loglog()
#~ plt.scatter(x_compare, y_compare, color='red')
#~ plt.plot(e_dist_x, e_dist_p)
#~ #plt.scatter(x, e_dist, facecolor='red', edgecolor='none')
plt.show()


