#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {


class HadronicInteraction: public Module {


public:
	HadronicInteraction( );
	void process(Candidate *candidate) const;
	double distribution_e(double energy, double x) const;
	double distribution_my1(double energy, double x) const; 
	double distribution_gamma(double energy, double x) const; 
	double distribution_Carceller(double energy, double x, double jcap, double a0, double b0) const;
	double distribution_Carceller_g(double energy, double x, double jcap, double a0, double b0) const;
	double CrossSection_Carceller(double energy) const;
	double CrossSection_Kelner(double energy) const;
	double CrossSection_Galprop(double energy) const;
	//~ double counter(Candidate *candidate, double density) const;
	//~ double counterPion(Candidate *candidate, double density) const;
};

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H
