from crpropa import *
from pylab import *
import random as rd
import numpy as np
import pandas as pd
import math



#Number of particles
particles = 1

#Steplength
steplength = 100000000000

#Candidate
c = Candidate()

h = HadronicInteraction()
#~ h.process(c)

F=h.distribution_gamma(1*TeV, 0.001)
print F
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
c.current.setEnergy(1 * TeV)
print c.current.getEnergy()
#~ m = ModuleList()
h.process(c)
#~ m.add(SimplePropagation(steplength * kpc,steplength* kpc))
h.process(c)
#~ m.add(HadronicInteraction())
#~ c.current.setEnergy(1 * TeV)
#~ c.setCurrentStep(10**5 * kpc)
#~ custep=c.getCurrentStep()/kpc
#~ energy=c.current.getEnergy()/TeV
#~ m.process(c)
#~ print custep
#~ print energy
h.process(c)
#~ m.add(obs)
#~ m.setShowProgress(True)
#~ m.run(source, particles, True)
