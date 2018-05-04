#include "crpropa/module/CreateElectrons.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"
#include "crpropa/module/MassDis.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {
    
    CreateElectrons::CreateElectrons() {
        
        setDescription("CreateElectrons");
    }
    
    void CreateElectrons::process(Candidate *candidate) const {
        double id = candidate->current.getId();
        double step = candidate->getTrajectoryLength();
        if (id=1000010010 and step < 1 * pc){
            crpropa::Vector3d pos= candidate->current.getPosition();
            candidate->addSecondary(11, 50*TeV, pos);
            
        }
        
    }
    
}
