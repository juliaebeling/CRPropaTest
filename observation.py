from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math

#Number of particles
particles = 500000000



#Candidate
c = Candidate()

#Radius of Sphere

r= 70* kpc



#Energy
E = 1000 * TeV
Energyvector=[]


#Steplength
steplength = 1.
source= Source()
source.add(SourceEnergy(E))
source.add(SourceParticleType(nucleusId(1,1)))
source.add(SourcePosition(Vector3d(0)))
output = TextOutput("Transfers/WithMassDis_Part"+str(particles)+"_SourceE"+str(E/TeV)+"_Step"+str(steplength)+"kpc_Obs"+str(r/ Mpc)+"Mpc_Sprotons.txt")
output.disableAll()
#output.enable(output.SerialNumberColumn)
###Setup turbulent magnetic field
randomSeed = 27
lMin=0.02 * pc
lMax=5*pc
l= turbulentCorrelationsLength(lMin, lMax, -11./3.)
spacing=0.025*pc/1.2
vgrid = VectorGrid(Vector3d(0), 1200, spacing)
vgrid.setReflective(True)
b= 350*1e-6*gauss
initTurbulence(vgrid, b, lMin, lMax, -11./3., randomSeed)
TurbField = MagneticFieldGrid(vgrid)
Magfiela.addField(TurbField)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.setEnergyScale(TeV)
obs=Observer()
obs.add(ObserverLargeSphere(Vector3d(0), r))
obs.onDetection(output)
m = ModuleList()
mod_PropCK=PropagationCK(MagField, 1e-3, steplength * kpc, steplength * kpc)
m.add(mod_PropCK)
m.add(HadronicInteraction())
m.add(obs)
m.add(MaximumTrajectoryLength(50*kpc))
m.setShowProgress(True)
m.run(source, particles, True)

