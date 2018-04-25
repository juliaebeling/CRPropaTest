from crpropa import *
from pylab import *
import random as rd
import numpy as np
import math

#Number of particles
particles = 10000



#Candidate
c = Candidate()

#Radius of Sphere

r= 250* pc




#Steplength
steplength = 0.1
source= Source()
origin = Vector3d(0)
for i in range(2600):
    source.add(SourceUniformCylinder(origin, 60*pc, 150*pc))
    source.add(SourceParticleType(nucleusId(1,1)))
    source.add(SourcePowerLawSpectrum(1*TeV, 1000*TeV, 2.0)
output = TextOutput("/Volumes/Samsung_T5/Studium/CRPropaTest/Transfers/NGC253.txt")
output.disableAll()
#output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.setEnergyScale(TeV)
obs=Observer()
obs.add(ObserverLargeSphere(Vector3d(0), r))
obs.onDetection(output)
               
    
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
m = ModuleList()
mod_PropCK=PropagationCK(MagField, steplength * pc, 10 * steplength * pc)
        #m.add(SimplePropagation(steplength*kpc, steplength*kpc))
m.add(HadronicInteraction())
m.add(EMInverseComptonScattering())
m.add(SynchrotronRadiation())
m.add(obs)
m.add(MaximumTrajectoryLength(1*kpc))
m.setShowProgress(True)
m.run(source, particles, True)


