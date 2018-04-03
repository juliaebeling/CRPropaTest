from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math

#Number of particles
particles = 1



#Candidate
c = Candidate()

#Radius of Sphere

r= 10* Mpc



#Energy
E = 10 * TeV
Energyvector=[]


#Steplength
steplength = 1000.
source= Source()
source.add(SourceEnergy(E))
source.add(SourceParticleType(nucleusId(1,1)))
source.add(SourcePosition(Vector3d(0)))
#output = TextOutput("/home/home1/jce/AnalyseCRPropa/TestGamma3_Part"+str(particles)+"_SourceE"+str(E/TeV)+"_Step"+str(steplength)+"kpc_Obs"+str(r/ Mpc)+"Mpc_Sprotons.txt")
output = TextOutput("/home/home1/jce/AnalyseCRPropa/1NewTest.txt")
output.disableAll()
output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.setEnergyScale(TeV)
obs=Observer()
obs.add(ObserverLargeSphere(Vector3d(0), r))
obs.onDetection(output)
RegField=UniformMagneticField(Vector3d(0, 0, 0)* 1e-6* gauss)
MagField=MagneticFieldList()
MagField.addField(RegField)
m = ModuleList()
mod_PropCK=PropagationCK(MagField, 1e-3, steplength * kpc, steplength * kpc)
m.add(mod_PropCK)
m.add(HadronicInteraction())
m.add(obs)
m.add(MaximumTrajectoryLength(10*Mpc))
m.setShowProgress(True)
m.run(source, particles, True)

