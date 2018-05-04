#ifndef CRPROPA_CREATEELECTRONS_H
#define CRPROPA_CREATEELECTRONS_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {
    
    
    class CreateElectrons: public Module {
        
        
        public:
        CreateElectrons( );
        void process(Candidate *candidate) const;
    };
    
} // namespace crpropa

#endif // CRPROPA_CREATEELECTRONS_H
