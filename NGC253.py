from crpropa import *
from pylab import *
import random as rd
import numpy as np
import math

#Number of particles
particles = 500



#Candidate
c = Candidate()

#Radius of Sphere

r= 162* pc




#Steplength
steplength = 0.1

sl=SourceList()
origin = Vector3d(0)

for i in range (2600):
    s=Source()
    phi = 2*M_PI*rd.random()
    RandRadius = 150*pc*pow(rd.random(), 1. / 2.)
    a= Vector3d (np.cos(phi)*RandRadius, np.sin(phi)*RandRadius, (-0.5+rd.random())*120*pc)
    s.add(SourcePosition(a))
    s.add(SourceParticleType(nucleusId(1,1)))
    #s.add(SourceParticleType(nucleusId(11,11)))
    s.add(SourcePowerLawSpectrum(1*TeV, 10000*TeV, 2.0))
    sl.add(s)





output = TextOutput("Simulations250418/TestEnergyChange.txt")
output.disableAll()
output.enable(output.SerialNumberColumn)
output.enable(output.CurrentIdColumn)
output.enable(output.CurrentEnergyColumn)
output.enable(output.CurrentPositionColumn)
output.enable(output.SourcePositionColumn)
output.setEnergyScale(TeV)
output.setLengthScale(pc)
obs=Observer()
obs.add(ObserverLargeSphere(Vector3d(0), r))
#obs.add(ObserverTimeEvolution(0*pc, 500*pc, 2))
#obs.setDeactivateOnDetection(False)
obs.onDetection(output)

MagField = MagneticFieldList()
randomSeed = 27
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
m.add(EMInverseComptonScattering(CMB, True, 1))
#.add(EMDoublePairProduction(CMB, True, 1))
m.add(EMPairProduction(CMB, True, 1))
#m.add(EMTripletPairProduction(CMB, True, 1))
#m.add(ElasticScattering(CMB))
m.add(ElectronPairProduction(CMB, True, 1))
#m.add(NuclearDecay())
#m.add(CreateElectrons())
m.add(SynchrotronRadiation(TurbField, True, 1))
#m.add(PhotoDisintegration(CMB, True, 1))
#m.add(PhotoPionProduction(CMB, True, True, False, 1, False))
m.add(obs)
m.add(MaximumTrajectoryLength(1*kpc))
m.setShowProgress(True)
m.run(sl, particles, True)


