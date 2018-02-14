#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {

HadronicInteraction::HadronicInteraction() {

	setDescription("HadronicInteraction");
}

double HadronicInteraction::distribution_e(double energy, double x) const{
	
	double L=log(energy / TeV);
	double Be= 1/(69.5+2.65*L+0.3*pow(L,2.));
	double betae=1/pow((0.201+0.062*L+0.00042*pow(L,2.)), 0.25);
	double ke=(0.279 + 0.141 *L + 0.0172* pow(L, 2.))/(0.3+ pow((2.3+L), 2.));
	double F=Be*pow((1+ke*pow(log(x),2.)), 3.) /(x*(1+0.3/pow(x, betae)))*(pow(-log(x), 5.));
	
	return F;
}

double HadronicInteraction::distribution_my1(double energy, double x) const{
	double L=log(energy / TeV);
	double Bm= 1.75+0.204*L+0.01 * pow(L,2.);
	double betam=1/(1.67+0.111*L+0.0038*pow(L,2.));
	double km=1.07-0.086*L+0.002*pow(L,2.);
	double aa=(1-pow(x,betam))/(1+km*pow(x, betam)*(1-pow(x,betam)));
	double A=Bm*log(x)/x*pow(aa, 4.);
	double B=1/log(x)-4*betam*pow(x,betam)/(1- pow(x,betam))-4*km*betam*pow(x, betam)*(1-2*pow(x,betam))/(1+km*pow(x,betam)*(1-pow(x,betam)));
	double F=A*B;
	
	return F;
}

double HadronicInteraction::distribution_gamma(double energy, double x) const{
	double L=log(energy / TeV);
	double Bg=1.3+0.14*L+0.011*L*L;
	double betag=1/(1.79+0.11*L+0.008*L*L);
	double kg=1/(0.801+0.049 *L+0.014 *L*L);
	double A=Bg*log(x)/x;
	double B=(1-pow(x, betag))/(1+kg*pow(x, betag)*(1-pow(x, betag)));
	double C=1/log(x)-4*betag*pow(x, betag)/(1-pow(x, betag))-4*kg*betag*pow(x, betag)*(1-2*pow(x, betag))/(1+kg*pow(x, betag)*(1-pow(x, betag)));
	double F=A*pow(B, 4.)*C;

	return F;
}

void HadronicInteraction::process(Candidate *candidate) const {
  double step = candidate->getCurrentStep();
  double energy = candidate->current.getEnergy();
  double id = candidate->current.getId();

  crpropa::Vector3d pos= candidate->current.getPosition();
  double cs_inel;
  double eG=energy / GeV;

    //~ calculates cross section for protons (energy dependend). Calculations based on Tan & Ng 1983
  if (id == 1000010010) {
	  double U = log(energy/ GeV * 1/200);
	 if (U >= 0 and energy >= 3 * GeV){
		 cs_inel=(32.2 * (1+0.0273*U))*1e-31+32.2*0.01*pow(U,2.)*1e-31;
	 }
	 if (U < 0 and energy >= 3 * GeV){
		cs_inel=(32.2 * (1+0.0273*U))*1e-31;
		 }
	 if (energy <= 0.3 * GeV){
		 cs_inel = 0;
	 }
	          }
	// Cross section based on Kelner 2009: 
	if (id == 1000010010){
		double L=log(energy / TeV);
		cs_inel=(34.3 + 1.88*L+0.25 *L*L)*1e-31;
	}
  //~ Cross section for different target material. Calculations based on R. Silberberg and C. H. Tsao
  //~ else {
	  //~ cs_inel=45 * pow(id, .7)*(1+0.016 * sin(5.3-2.63*log(id)))*1e-31 ;
	  //~ cs_inel=cs_inel*(1-0.62* pow(exp(1), -energy/200)*sin(10.9 * pow(10.9*energy, -0.28)));
		//~ }


Random &random = Random::instance();
double p_pp=cs_inel*1e6*step;

double ra = random.rand();


  if (ra > p_pp){

	return;
	  }

double limit = 1 / p_pp;

//~ if (step > limit) {
			//~ // limit next step to mean free path
			//~ candidate->limitNextStep(limit);
	
//~ }
double Eout=0;
label:
double test=1;
double Emo;
double Ee;
double Ene;
double Emt;
double Eg;
double gamma=0;
//~ Gamma rays

	do{
	double x=random.rand()*(-3);
	double F=distribution_gamma(energy, pow(10, x));
	double Fmax=distribution_gamma(energy, 0.001);
	double y=random.rand()*Fmax;
		if (y < F){
			double Eout=pow(10, x)*energy;
			candidate->addSecondary(22, Eout, pos);
			gamma++;
			Eg=Eout;
			}
		}while (gamma ==0);


//~ First myon neutrino
	do {
	double x=random.rand()*(-3);
	double F=distribution_my1(energy, pow(10, x));
	double Fmax=distribution_my1(energy, 0.001);
	double y=random.rand()*Fmax;
		if (y < F){
			double Eout=pow(10, x)*energy;
			candidate->addSecondary(14, Eout, pos);
			test=2;
			Emo=Eout;
		}
		} while (test == 1);

//~ Electron
	do {
	double x=random.rand()*(-3);

	double F=distribution_e(energy, pow(10, x));
	double k=0.001;
	double Fmax=distribution_e(energy, 0.001);
	double Eout=F;
	double y=random.rand()*Fmax;
	Eout=pow(10, x)*energy;
		if (y < F and (Eout+Emo)<energy){
			candidate->addSecondary(11, Eout, pos);
			test=3;
			Ee=Eout;

					}
		} while (test == 2);

//~ Electron neutrino
	do {
	double x=random.rand()*(-3);

	double F=distribution_e(energy, pow(10, x));
	double Fmax=distribution_e(energy, 0.001);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F and (Eout+Emo)<energy){
			candidate->addSecondary(12, Eout, pos);
			test=4;
			Ene=Eout;
				}

		} while (test == 3);


//~ Second myon neutrino
		do {
	double x=random.rand()*(-3);


	double F=distribution_e(energy, pow(10, x));
	double Fmax=distribution_e(energy, 0.001);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F and (Eout+Emo)<energy){
			candidate->addSecondary(14, Eout, pos);
			test=5;
					}
		} while (test == 4);

	if (Ee+Ene+Emt+Emo+Eg > energy)
		{
		std::cout << "H";
		goto label;
		}
		return;
}




} //~ namespace CRPropa























//~ let's particle propagate until it interacts. returns traveled distance until interaction

//~ double HadronicInteraction::counter(Candidate *candidate, double density) const {
	
  //~ double step = candidate->getCurrentStep();
  //~ double energy = candidate->current.getEnergy();
  //~ initializes return value (count)
  //~ double count = 0;
          //~ particle propagates until it has interacted
          //~ do {  
		  //~ double step = candidate->getCurrentStep();
	      //~ double U = log(energy/ GeV * 1/200);
	      //~ double cs_inel=(32.2 * (1+0.0273*U))*1e-31;
	         //~ if ( U >= 0.0){
	                 //~ cs_inel=cs_inel+32.2*0.01*pow(U,2.)*1e-31;
	                       //~ }
		//~ Random &random = Random::instance();
		//~ double p_pp=cs_inel*density*step;
		//~ double ra = random.rand();
		  //~ if (ra < p_pp){
		  //~ candidate->setActive(false);
		  //~ }
		  //~ count = count + step;
	     //~ }while (candidate->isActive() == true);


	//~ double co= count /kpc;
	//~ return co;
//~ }

//~ double HadronicInteraction::counterPion(Candidate *candidate, double density) const {
	
  //~ double step = candidate->getCurrentStep();
  //~ double energy = candidate->current.getEnergy();
  //~ crpropa::Vector3d pos= candidate->current.getPosition();
  //~ initializes return value (count)
  //~ double count = 0;
          //~ particle propagates until it has interacted
          //~ do {  
		 //~ double p_nnpp=9e-32*density*step;
		 //~ double p_pppm=2e-31*density*step;
		 //~ double p_ppnn=3e-32*density*step;
		 //~ double p_pnpn=2e-31*density*step;
		 //~ double Eout = 1 *  TeV; 
		 //~ int sign = 1;
		//~ Random &random = Random::instance();
		//~ double ra = random.rand();
		  //~ if (ra < p_nnpp){
		  //~ candidate->setActive(false);
		  //~ //nu_myon
		  //~ candidate->addSecondary(sign * 14, Eout, pos);
		  //~ candidate->addSecondary(sign * 14, Eout, pos);
		  //~ //antinu_myon
		  //~ candidate->addSecondary(sign * -14, Eout, pos);
		  //~ candidate->addSecondary(sign * -14, Eout, pos);
		  //~ //nu_e
		  //~ candidate->addSecondary(sign * 12, Eout, pos);
		  //~ candidate->addSecondary(sign * 12, Eout, pos);
		  //~ //electron
		  //~ candidate->addSecondary(sign * 11, Eout, pos);
		  //~ candidate->addSecondary(sign * 11, Eout, pos);
		  //~ }
		  //~ if ( p_nnpp <= ra and ra < (p_nnpp+p_pppm)){
		  //~ candidate->setActive(false);
		  //~ //electron
		  //~ candidate->addSecondary(sign * 11, Eout, pos);
		  //~ //antinu_e
		  //~ candidate->addSecondary(sign * -12, Eout, pos);
		  //~ //nu_myon
		  //~ candidate->addSecondary(sign * 14, Eout, pos);
		  //~ //nu_e
		  //~ candidate->addSecondary(sign * 12, Eout, pos);
		  //~ //antinu_myon
		  //~ candidate->addSecondary(sign * -14, Eout, pos);
		  //~ //positron
		  //~ candidate->addSecondary(sign * -11, Eout, pos);
		  //~ }
		  //~ if ((p_nnpp+p_pppm) <= ra and ra < (p_nnpp+p_pppm+p_ppnn)){
		  //~ candidate->setActive(false);
		  //~ //photon
		  //~ candidate->addSecondary(22, Eout, pos);
		  //~ candidate->addSecondary(22, Eout, pos);
		  //~ }
		  //~ if ((p_nnpp+p_pppm+p_ppnn) <= ra and ra < (p_nnpp+p_pppm+p_ppnn+p_pnpn)){
		  //~ candidate->setActive(false);
		  //~ //positron
		  //~ candidate->addSecondary(sign * -11, Eout, pos);
		  //~ //nu_e
		  //~ candidate->addSecondary(sign * 12, Eout, pos);
		  //~ //antinu_myon
		  //~ candidate->addSecondary(sign * -14, Eout, pos);
		  //~ //nu_myon
		  //~ candidate->addSecondary(sign * 14, Eout, pos);
		  //~ //photon
		  //~ candidate->addSecondary(22, Eout, pos);
		  //~ }
		  //~ count = count + step;
	     //~ }while (candidate->isActive() == true);

	//~ double co= count /kpc;
	//~ return co;
//~ }

