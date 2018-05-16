from crpropa import *
from pylab import *
import random as rd
import numpy as np
import math
HI=HadronicInteraction()
#Number of particles
particles = 10000


#Candidate
c = Candidate()

#Radius of Sphere

r= 200*pc

#Steplength
steplength = 1

#Initialise Sourcelist
sl=SourceList()
origin = Vector3d(0)

#Create Sources
#Calculated by Integeral
for i in range (50):
    s=Source()
    s.add(SourcePosition(HI.Position(120*pc,200*pc)))
    s.add(SourceParticleType(nucleusId(1,1)))
    s.add(SourcePowerLawSpectrum(1*TeV, 10000*TeV, -2.0))
    sl.add(s)

#Create output1
output = TextOutput("Simulations150518/1_OS.txt")
output.disableAll()
output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.enable(output.CurrentPositionColumn)
output.enable(output.SourcePositionColumn)
output.setEnergyScale(TeV)
output.setLengthScale(pc)

#Create Observer Sphere
obs=Observer()
obs.add(ObserverLargeSphere(origin, r))
obs.onDetection(output)


#Create output2
output2 = TextOutput("Simulations150518/1_Source.txt")
output2.disableAll()
output2.enable(output.SerialNumberColumn)
output2.enable(output.CurrentIdColumn)
output2.enable(output.CurrentEnergyColumn)
output2.enable(output.CurrentPositionColumn)
output2.enable(output.SourcePositionColumn)
output2.setEnergyScale(TeV)
output2.setLengthScale(pc)

#Create Time Evolution Source
obs2=Observer()
obs2.add(ObserverTimeEvolution(0, 10*kpc, 1))
obs2.setDeactivateOnDetection(False)
obs2.onDetection(output2)

#Create output3 5 kpc
output3 = TextOutput("Simulations150518/1_5kpc.txt")
output3.disableAll()
output3.enable(output.SerialNumberColumn)
output3.enable(output.CurrentIdColumn)
output3.enable(output.CurrentEnergyColumn)
output3.enable(output.CurrentPositionColumn)
output3.enable(output.SourcePositionColumn)
output3.setEnergyScale(TeV)
output3.setLengthScale(pc)

#Create Time Evolution kpc
obs3=Observer()
obs3.add(ObserverTimeEvolution(5*kpc, 10*kpc, 1))
obs3.setDeactivateOnDetection(False)
obs3.onDetection(output3)


#Initiate turbulent magnetic field
MagField = MagneticFieldList()
randomSeed = 23
lMin=0.04 * pc
lMax=2.*pc
l= turbulentCorrelationLength(lMin, lMax, -11./3.)
spacing=0.01*pc
vgrid = VectorGrid(origin, 600, spacing)
vgrid.setReflective(True)
b= 10*1e-6*gauss
initTurbulence(vgrid, b, lMin, lMax, -11./3., randomSeed)
TurbField = MagneticFieldGrid(vgrid)
MagField.addField(TurbField)

#ModuleList
m = ModuleList()
#Propagation
mod_PropCK=PropagationCK(MagField, 0.5, steplength * pc, 100 * steplength * pc)
m.add(mod_PropCK)

#Interaction
m.add(HadronicInteraction())
#m.add(EMInverseComptonScattering(CMB, True, 1))
#m.add(EMDoublePairProduction(CMB, True, 1))
#m.add(EMPairProduction(CMB, True, 1))
#m.add(EMTripletPairProduction(CMB, True, 1))
#m.add(ElasticScattering(CMB))
#m.add(ElectronPairProduction(CMB, True, 1))
#m.add(SynchrotronRadiation(TurbField, True, 1))
#m.add(PhotoDisintegration(CMB, True, 1))
#m.add(PhotoPionProduction(CMB, True, True, False, 1, False))

#Observer
m.add(obs)
m.add(obs2)
m.add(obs3)
m.add(MaximumTrajectoryLength(6*kpc))

m.setShowProgress(True)
m.run(sl, particles, True)





