#ifndef CRPROPA_MASSDIS_H
#define CRPROPA_MASSDIS_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {
    
    
    class MassDis: public Module {
        
        
        public:
        MassDis( );
        static double getDensity(double x, double y, double z);
        static double NGC253(double x, double y, double z);
    };
    
} // namespace crpropa

#endif // CRPROPA_MASSDIS_H
