#include "crpropa/module/MassDis.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {
    MassDis::MassDis() {
        
        setDescription("MassDis");
    }

double MassDis::getDensity(double x, double y, double z) {
    
    double R = sqrt(pow(x,2)+pow(y,2));
    double ht=0.15*kpc*exp((R-9.8*kpc)/R);
    
    
    if(R<7*kpc) //INNEN
    {
        double density= 1.672621898*pow(10, -27)*1.5*pow(10,6)*exp(-(7*kpc-8.4*kpc)/(3.15*kpc)-z*log(2)/ht);
        return density;
    }
    else{ //AUßEN
        double density= 1.672621898*pow(10, -27)*1.5*pow(10,6)*exp(-(R-8.4*kpc)/(3.15*kpc)-z*log(2)/ht);
        return density;
        
        }
    
}

    double MassDis::NGC253(double x, double y, double z){
        double density=0;
        if(std::abs(x) < 150*pc and std::abs(y) < 150*pc and std::abs(y) < 60 *pc)
        {
            density= 1160*pow(10, 6);
        }
        return density;
    }
    
}//END NAMESPACE CRPROPA
