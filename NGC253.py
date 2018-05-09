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

r= 162*pc




#Steplength
steplength = 0.1

sl=SourceList()
origin = Vector3d(0)
b= Vector3d(0)
for i in range (2600):
    s=Source()
    phi = 2*M_PI*rd.random()
    RandRadius = 150*pc*pow(rd.random(), 1. / 2.)
    a= Vector3d (np.cos(phi)*RandRadius, np.sin(phi)*RandRadius, (-0.5+rd.random())*120*pc)
    s.add(SourcePosition(a))
    s.add(SourceParticleType(nucleusId(1,1)))
    s.add(SourcePowerLawSpectrum(1*TeV, 10000*TeV, 2.0))
    sl.add(s)





output = TextOutput("Simulations070518/onlyHIsource.txt")
output2 = TextOutput("Simulations070518/onlyHIObserverSphere.txt")
output.disableAll()
output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.enable(output.CurrentPositionColumn)
output.enable(output.SourcePositionColumn)
output.setEnergyScale(TeV)
output.setLengthScale(pc)
output2.disableAll()
output2.enable(output.SerialNumberColumn)
output2.enable(output.CurrentIdColumn)
output2.enable(output.CurrentEnergyColumn)
output2.enable(output.CurrentPositionColumn)
output2.enable(output.SourcePositionColumn)
output2.setEnergyScale(TeV)
output2.setLengthScale(pc)
obs=Observer()
obs2=Observer()
obs.add(ObserverLargeSphere(b, r))
#obs.add(ObserverTimeEvolution(10*kpc, 1*kpc, 1))
obs2.add(ObserverTimeEvolution(0*kpc, 5*kpc, 1))
obs.setDeactivateOnDetection(False)
obs2.setDeactivateOnDetection(False)
obs.onDetection(output)
obs2.onDetection(output2)

MagField = MagneticFieldList()
randomSeed = 23
lMin=0.04 * pc
lMax=2.*pc
l= turbulentCorrelationLength(lMin, lMax, -11./3.)
spacing=0.01*pc
vgrid = VectorGrid(Vector3d(0), 600, spacing)
vgrid.setReflective(True)
b= 350*1e-6*gauss
initTurbulence(vgrid, b, lMin, lMax, -11./3., randomSeed)
TurbField = MagneticFieldGrid(vgrid)
MagField.addField(TurbField)
m = ModuleList()
mod_PropCK=PropagationCK(MagField, 0.5, steplength * pc, 100 * steplength * pc)
#m.add(SimplePropagation(steplength*kpc, steplength*kpc))
m.add(mod_PropCK)
m.add(HadronicInteraction())
#m.add(EMInverseComptonScattering(CMB, True, 1))
#m.add(EMDoublePairProduction(CMB, True, 1))
#m.add(EMPairProduction(CMB, True, 1))
#m.add(EMTripletPairProduction(CMB, True, 1))
#m.add(ElasticScattering(CMB))
#m.add(ElectronPairProduction(CMB, True, 1))
#m.add(NuclearDecay())
#m.add(CreateElectrons())
#m.add(SynchrotronRadiation(TurbField, True, 1))
#m.add(PhotoDisintegration(CMB, True, 1))
#m.add(PhotoPionProduction(CMB, True, True, False, 1, False))
m.add(obs)
m.add(obs2)
m.add(MaximumTrajectoryLength(4.0001*kpc))
m.setShowProgress(True)
m.run(sl, particles, True)


