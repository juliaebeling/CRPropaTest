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

//Energy distribution for pions

double HadronicInteraction::distribution_pi(double energy, double x) const{
	double L=log(energy/TeV);
	double Bp=5.58+0.78*L+0.1*L*L;
	double r=3.1/pow(Bp, 1.5);
	double a=0.89/(pow(Bp, 0.5)*(1-exp(-0.33*Bp)));
	double mp=134.9766*MeV;
	double A=4*a*Bp*pow(x, a-1);
	double B=(1-pow(x,a))/(1+r*pow(x,a)*(1-pow(x,a)));
	double C=(1/(1-pow(x,a))+(r*(1-2*pow(x,a))/(1+r*pow(x,a)*(1-pow(x,a)))))*pow(1-mp/(x*energy), 0.5);
	double F=A*pow(B, 4.)*C;
	
	return F;
}

double HadronicInteraction::number_pi(double energy) const{
	double x=1/1000.;
	double i=1/100000.;
	double y=0;
	double j=0;
	do{
	y=y+distribution_pi(energy,x);
	x=x+i;
	j++;
	}while(x < 1);
	return y/j*(x-1/1000.);
}

//Energy distribution for electron, electron neutrino and second muon neutrino based on Kelner 2006

double HadronicInteraction::distribution_e(double energy, double x) const{
	
	double L=log(energy / TeV);
	double Be= 1/(69.5+2.65*L+0.3*pow(L,2.));
	double betae=1/pow((0.201+0.062*L+0.00042*pow(L,2.)), 0.25);
	double ke=(0.279 + 0.141 *L + 0.0172* pow(L, 2.))/(0.3+ pow((2.3+L), 2.));
	double F=Be*pow((1+ke*pow(log(x),2.)), 3.) /(x*(1+0.3/pow(x, betae)))*(pow(-log(x), 5.));
	
	return F;
}

double HadronicInteraction::number_e(double energy) const{
	double x=1/1000.;
	double i=1/100000.;
	double y=0;
	double j=0;
	do{
	y=y+distribution_e(energy,x);
	x=x+i;
	j++;
	}while(x < 1);
	return y/j*(x-1/1000.);
}



//Energy distribution for first muon neutrino based on Kelner 2006
double HadronicInteraction::distribution_my1(double energy, double x) const{
	double L=log(energy / TeV);
	double Bm= 1.75+0.204*L+0.01 * pow(L,2.);
	double betam=1/(1.67+0.111*L+0.0038*pow(L,2.));
	double km=1.07-0.086*L+0.002*pow(L,2.);
	x=x/0.427;
	double aa=(1-pow(x,betam))/(1+km*pow(x, betam)*(1-pow(x,betam)));
	double A=Bm*log(x)/x*pow(aa, 4.);
	double B=1/log(x)-4*betam*pow(x,betam)/(1- pow(x,betam))-4*km*betam*pow(x, betam)*(1-2*pow(x,betam))/(1+km*pow(x,betam)*(1-pow(x,betam)));
	double F=A*B;
	
	return F;
}

double HadronicInteraction::number_my1(double energy) const{
	double x=1/1000.;
	double i=1/100000.;
	double y=0;
	double j=0;
	do{
	y=y+distribution_my1(energy,x);
	x=x+i;
	j++;
	}while(x < 1);
	return y/j*(x-1/1000.);
}
//Energy distribution for gamma photons based on Kelner 2006
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

double HadronicInteraction::number_gamma(double energy) const{
	double x=1/1000.;
	double i=1/100000.;
	double y=0;
	double j=0;
	do{
	y=y+distribution_gamma(energy,x);
	x=x+i;
	j++;
	}while(x < 1);
	return y/j*(x-1/1000.);
}
//Energy distribution for lepton secondaries of pp interactions based on Carceller 2017

double HadronicInteraction::distribution_Carceller(double energy, double x, double jcap, double a0, double b0) const{
	double a = a0 * (1+ 0.073*log(energy/ PeV)+ 0.0070*log(energy/PeV)*log(energy/PeV)); 
	double b = b0 * (1+0.020*log(energy/PeV)+0.0018*log(energy/PeV)*log(energy/PeV));
	double A = a * pow((1-jcap*x), 3.)/x;
	double B = exp(-b*pow(jcap*x, 0.43))/pow(1+pow(0.1*GeV/(x*energy), 0.5), 2.);
	double F=A*B; 
	
	return F; 
}


//Energy distribution for gamma photons based on Carceller 2017

double HadronicInteraction::distribution_Carceller_g(double energy, double x, double jcap, double a0, double b0) const{
	double a = a0 * (1+ 0.073*log(energy/ PeV)+ 0.0070*log(energy/PeV)*log(energy/PeV)); 
	double b = b0 * (1+0.02*log(energy/PeV)+0.0018*log(energy/PeV)*log(energy/PeV));
	double A = a * pow((1-jcap*x), 3.)/x;
	double B = exp(-b*pow(jcap*x, 0.43))/pow(1+pow(0.2*GeV/(x*energy), 0.5), 2.);
	double F=A*B; 
	
	return F; 
}

//Cross Section of inelastic pp interaction based on Tan & Ng 1983 (Used in Galprop)
double HadronicInteraction::CrossSection_Galprop(double energy) const{
	double cs_inel;
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
	 return cs_inel;
	}

//Cross Section of inelastic pp interaction based on Kelner 2006
double HadronicInteraction::CrossSection_Kelner(double energy) const{
		double L=log(energy / TeV);
		double A=1-pow(1.22*1e-3*TeV/energy, 4.);
		double cs_inel=(34.3 + 1.88*L+0.25 *L*L)*A*A*1e-31;
		return cs_inel;
}

//Cross Section of inelastic pp interaction based on Carceller 2017
double HadronicInteraction::CrossSection_Carceller(double energy) const{
	double cs_inel=17.7*pow(energy/GeV, 0.082)*1e-31;
	return cs_inel;
}

//Process Function:

void HadronicInteraction::process(Candidate *candidate) const {

  double step = candidate->getCurrentStep();
  double energy = candidate->current.getEnergy();
  double id = candidate->current.getId();
  //~ if (energy < 10 * TeV){
	  //~ return;
	  //~ }

if (id < 500){
	return;
			}

  crpropa::Vector3d pos= candidate->current.getPosition();
  double cs_inel=0;
  double eG=energy / GeV;
  double jcap=1;


if (id == 1000010010) {
	  cs_inel=CrossSection_Galprop(energy);
					  }

if (id == 1000010010) {
	cs_inel=CrossSection_Kelner(energy);
					  }


	//~ if (id == 1000010010){
	//~ cs_inel=CrossSection_Carceller(energy);
						 //~ }


if (id == 1000020040 and energy > 1* TeV) {
	double cs_Hep=60.5*1e-31;
	cs_inel= cs_Hep*pow(eG, 0.062);
	jcap=4;
	}


if (id == 1000260056 and energy > 1* TeV) {
	double cs_Fep=551*1e-31;
	cs_inel= cs_Fep*pow(eG, 0.026);
	jcap=56; 
	}

  //~ Cross section for different target material. Calculations based on R. Silberberg and C. H. Tsao
  //~ else {
	  //~ cs_inel=45 * pow(id, .7)*(1+0.016 * sin(5.3-2.63*log(id)))*1e-31 ;
	  //~ cs_inel=cs_inel*(1-0.62* pow(exp(1), -energy/200)*sin(10.9 * pow(10.9*energy, -0.28)));
		//~ }


Random &random = Random::instance();

double density=1e6;
double p_pp=cs_inel*density*step;
double ra = random.rand();

  if (ra > p_pp or energy < 1*GeV){

	return;

	  }


double limit = 1 / p_pp*0.1;

if (step > limit) {
			// limit next step to mean free path
			candidate->limitNextStep(limit);
	
}
double Eout=0;


double Emo;
double Ee;
double Ene;
double Emt;
double Eg;
double gamma=0;
double test=0;
//~ goto label2;
if (jcap == 1)
{
	goto label; 
}


label:

//~ Gamma rays


	
	double r=random.rand();
	double F=number_gamma(energy);
	std::cout << F << std::endl; 
	double N=std::floor(F);
	double Econ=0;
	
	if (r < F-N){
		N=N+1;
	}

	double i=1;
		do{
			double x=random.rand()*(-3);
			 Eout=pow(10, x)*energy;
			double E=distribution_gamma(energy, pow(10, x));
			double Emax=distribution_gamma(energy, 0.001);
			double y=random.rand()*Emax;
			double Eoutg=pow(10, x)*energy;
			double Et=Econ+Eout;
			if (y < E and Et < energy){
				candidate->addSecondary(22, Eout, pos);
				Econ=Econ+Eout;
				i++;
				std::cout << Econ/energy << std::endl;}
			}while (i <= N);
		
	Eg=Econ;
	Econ=0;



//~ First myon neutrino 14

	r=random.rand();
	F=number_my1(energy);
	std::cout << F << std::endl;
	N=std::floor(F);
	if (r < F-N){
		N=N+1;
	}
	i=1;
			do{
			double x=random.rand()*(-3);
			 Eout=pow(10, x)*energy;
			double E=distribution_my1(energy, pow(10, x));
			double Emax=distribution_my1(energy, 0.001);
			double y=random.rand()*Emax;
			double Et=Econ+Eout;
			double Eout=pow(10, x)*energy;
			if (y < E and Et< (energy-Eg)){
				candidate->addSecondary(14, Eout, pos);
				Econ=Econ+Eout,
				std::cout << Econ/energy << std::endl;
				i++;}
			}while (i <= N);
	Emo=Econ;
	Econ=0;


//~ Electron 22

	r=random.rand();
	F=number_e(energy);
	std::cout << F << std::endl;
	N=std::floor(F);
	if (r < F-N){
		N=N+1;
	}
	i=1;
		do{
		double x=random.rand()*(-3);
		 Eout=pow(10, x)*energy;
		double E=distribution_e(energy, pow(10, x));
		double Emax=distribution_e(energy, 0.001);
		double y=random.rand()*Emax;
		double Et=Econ+Eout;
		double Eout=pow(10, x)*energy;
		if (y < E and Et < (energy-Emo-Eg)){
			candidate->addSecondary(22, Eout, pos);
			Econ=Econ+Eout;
			std::cout << Econ/energy << std::endl;
			i++;}
		}while (i <= N);
	Ee=Econ;
	Econ=0;


//~ Electron neutrino 12

	r=random.rand();
	std::cout << F << std::endl;

	if (r < F-N){
		N=N+1;
	}
	i=1;
	do{
		double x=random.rand()*(-3);
		 Eout=pow(10, x)*energy;
		double E=distribution_e(energy, pow(10, x));
		double Emax=distribution_e(energy, 0.001);
		double y=random.rand()*Emax;
		double Et=Econ+Eout;
		double Eout=pow(10, x)*energy;
		if (y < E and Et < (energy-Emo-Eg-Ee)){
			candidate->addSecondary(12, Eout, pos);
			Econ=Econ+Eout;
			std::cout << Econ/energy << std::endl;
			i++;}
		}while (i <= N);
	Ene=Econ;
	Econ=0;


//~ Second myon neutrino 14

	r=random.rand();
	std::cout << F << std::endl;

	if (r < F-N){
		N=N+1;
	}
	i=1;
	do{
		double x=random.rand()*(-3);
		 Eout=pow(10, x)*energy;
		double E=distribution_e(energy, pow(10, x));
		double Emax=distribution_e(energy, 0.001);
		double y=random.rand()*Emax;
		double Eout=pow(10, x)*energy;
		if (y < F and Econ < (energy-Emo-Eg-Ee-Ene)){
			candidate->addSecondary(14, Eout, pos);
			Econ=Econ+Eout;
			std::cout << 's' << Econ/energy << std::endl;
			i++;}
		}while (i <= N);
	Emt=Econ;
	Econ=0;
		
	if (Ee+Ene+Emt+Emo+Eg> energy)
		{
		std::cout << "H";

		}
	candidate->current.setEnergy(energy-(Ee+Ene+Emt+Emo+Eg));
	std::cout << (Ee+Ene+Emt+Emo+Eg)/energy << std::endl;
	return;


label2: 
test = 1; 

//gamma
	do {
	double x=random.rand()*(-3);


	double F=distribution_Carceller_g(energy, pow(10, x), jcap, 5.8, 5.3);
	double Fmax=distribution_Carceller_g(energy, 0.001, jcap, 5.8, 5.3);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F ){
			candidate->addSecondary(22, Eout, pos);
			test=2;
					}
		} while (test == 1);

//electron neutrino

	do {
	double x=random.rand()*(-3);


	double F=distribution_Carceller(energy, pow(10, x), jcap, 2.7, 7.7);
	double Fmax=distribution_Carceller(energy, 0.001, jcap, 2.7, 7.7);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F ){
			candidate->addSecondary(12, Eout, pos);
			test=3;
					}
		} while (test == 2);

//antielectron neutrino

	do {
	double x=random.rand()*(-3);


	double F=distribution_Carceller(energy, pow(10, x), jcap, 2.7, 8.7);
	double Fmax=distribution_Carceller(energy, 0.001, jcap, 2.7, 8.7);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F ){
			candidate->addSecondary(-12, Eout, pos);
			test=4;
					}
		} while (test == 3);

//muon neutrino

	do {
	double x=random.rand()*(-3);


	double F=distribution_Carceller(energy, pow(10, x), jcap, 5.1, 8.3);
	double Fmax=distribution_Carceller(energy, 0.001, jcap, 5.1, 8.3);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F ){
			candidate->addSecondary(14, Eout, pos);
			test=4;
					}
		} while (test == 3);

//antimuon neutrino
	do {
	double x=random.rand()*(-3);


	double F=distribution_Carceller(energy, pow(10, x), jcap, 5.1, 8.3);
	double Fmax=distribution_Carceller(energy, 0.001, jcap, 5.1, 8.3);
	double y=random.rand()*Fmax;
	double Eout=pow(10, x)*energy;
		if (y < F ){
			candidate->addSecondary(-14, Eout, pos);
			test=5;
					}
		} while (test == 4);
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

