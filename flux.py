from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math

e_dist_my=[]
e_dist_e=[]
e_dist_ne=[]
e_dist_g=[]

data1=[]
data1 = pd.read_csv('Part1000000SourceE1000.0MaxStep10protons.txt',
       names=['SN','ID','E','SNO','SN1'], delimiter='\t', comment='#',
       usecols=["SN", "ID", "E"])
#~ data.sort_values('E', ascending=True)
data1['ID']=data1.ID
data1['E']=data1.E

data1.ID.tolist()
data1.E.tolist()


data2=[]
data2 = pd.read_csv('Part10000000_SourceE1.0_MaxStep100kpc_Obs10.0Mpc_Sprotons.txt'
       names=['SN','ID','E','SNO','SN1'], delimiter='\t', comment='#',
       usecols=["SN", "ID", "E"])
#~ data.sort_values('E', ascending=True)
data2['ID']=data2.ID
data2['E']=data2.E

data2.ID.tolist()
data2.E.tolist()

e_dist2_my=[]
e_dist2_e=[]
e_dist2_ne=[]
e_dist2_g=[]
#Einlesen der verschiedenen Teilchen
for i in range (len(data1.ID)):
	#~ if data1.ID[i]==11:
		#~ e_dist2_e.append(data1.E[i]*EeV)
	#~ if data1.ID[i]==14 or data1.ID[i]==-14 :
		#~ e_dist2_my.append(data1.E[i]*EeV)
	#~ if data1.ID[i]==12 or data1.ID[i]==-12 :
		#~ e_dist2_ne.append(data1.E[i]*EeV)
	if data1.ID[i]==22:
		e_dist_g.append(data1.E[i]*EeV)
		

#Einlesen der verschiedenen Teilchen
for i in range (len(data2.ID)):
	#~ if data2.ID[i]==11:
		#~ e_dist_e.append(data2.E[i]*EeV)
	#~ if data2.ID[i]==14 or data2.ID[i]==-14 :
		#~ e_dist_my.append(data1.E[i]*EeV)
	#~ if data2.ID[i]==12 or data2.ID[i]==-12 :
		#~ e_dist_ne.append(data2.E[i]*EeV)
	if data2.ID[i]==22:
		e_dist2_g.append(data2.E[i]*EeV)


#Number of bins
nob=100

#Edges for histogram bins 1

Edges = np.logspace(-10, -3, nob)
h, binEdges_e = np.histogram(e_dist_e, bins=Edges)
my, binEdges_my = np.histogram(e_dist_my, bins=Edges)
g, binEdges_g=np.histogram(e_dist_g, bins=Edges)
binMiddle_e = (binEdges_e[1:]-binEdges_e[:-1])/2.+binEdges_e[:-1]
binMiddle_my = (binEdges_my[1:]-binEdges_my[:-1])/2.+binEdges_my[:-1]
binMiddle_g=(binEdges_g[1:]-binEdges_g[:-1])/2.+binEdges_g[:-1]

#Edges for histogram bins 2
Edges2 = np.logspace(-10, -3, nob)
h2, binEdges2_e = np.histogram(e_dist2_e, bins=Edges2)
my2, binEdges2_my = np.histogram(e_dist2_my, bins=Edges2)
g2, binEdges2_g=np.histogram(e_dist2_g, bins=Edges2)
binMiddle2_e = (binEdges2_e[1:]-binEdges2_e[:-1])/2.+binEdges2_e[:-1]
binMiddle2_my = (binEdges2_my[1:]-binEdges2_my[:-1])/2.+binEdges2_my[:-1]
binMiddle2_g=(binEdges2_g[1:]-binEdges2_g[:-1])/2.+binEdges2_g[:-1]



#~ plt.plot(binMiddle_e, h, linewidth=0., marker='s', label = 'Electrons')
#~ plt.plot(binMiddle_my, my,linewidth=0., marker='s', label='Muon Neutrinos')
plt.plot(binMiddle_g, g,linewidth=0., marker='s', label='Gamma rays')

#~ plt.plot(binMiddle2_e, h2, linewidth=0., marker='s', label = 'Electrons2')
#~ plt.plot(binMiddle2_my, my2,linewidth=0., marker='s', label='Muon Neutrinos2')
plt.plot(binMiddle2_g, g2, label='Gamma rays2')

plt.legend(loc='upper right')
plt.loglog()
#~ plt.scatter(x_compare, y_compare, color='red')
#~ plt.plot(e_dist_x, e_dist_p)
#~ #plt.scatter(x, e_dist, facecolor='red', edgecolor='none')
plt.show()
fig=plt.figure()
