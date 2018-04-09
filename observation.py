from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math

#Number of particles
particles = 100000



#Candidate
c = Candidate()

#Radius of Sphere

r= 10* Mpc



#Energy
E = 1000 * TeV
Energyvector=[]


#Steplength
steplength = 100.
source= Source()
source.add(SourceEnergy(E))
source.add(SourceParticleType(nucleusId(1,1)))
source.add(SourcePosition(Vector3d(0)))
output = TextOutput("Transfers/CC2ompareKelner_Part"+str(particles)+"_SourceE"+str(E/TeV)+"_Step"+str(steplength)+"kpc_Obs"+str(r/ Mpc)+"Mpc_Sprotons.txt")
output.disableAll()
#output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.setEnergyScale(TeV)
obs=Observer()
obs.add(ObserverLargeSphere(Vector3d(0), r))
obs.onDetection(output)
#RegField=UniformMagneticField(Vector3d(0, 0, 0)* 1e-6* gauss)
#MagField=MagneticFieldList()
#MagField.addField(RegField)
m = ModuleList()
#mod_PropCK=PropagationCK(MagField, 1e-3, steplength * kpc, steplength * kpc)
m.add(SimplePropagation(steplength*kpc, steplength*kpc))
m.add(HadronicInteraction())
m.add(obs)
m.add(MaximumTrajectoryLength(20*Mpc))
m.setShowProgress(True)
m.run(source, particles, True)

