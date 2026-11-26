#include "DoubleNModel.hpp"
#include <Tools/Core/Spinor.hpp>
#include <Glauber/MeanFieldSpinor.hpp>
#include <Tools/Core/Utilfunctions.hpp>

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
using namespace nucqcd;

using namespace std;



DoubleNModel::DoubleNModel(MeanFieldNucleusThick *pnucleus,bool setpw, 
			   bool setSRC, bool setCT, bool setcorr, int type1, int type2, 
			   double precision, int integr, string dir)
:pw(setpw), SRC(setSRC), CT(setCT), corr(setcorr), 
particletype1(type1), particletype2(type2), pnucl(pnucleus), prec(precision), integrator(integr), homedir(dir){
}


DoubleNModel::~DoubleNModel(){
  
}


complex<double> DoubleNModel::getMatrixEl(const Kinematics2to3 &tk, int spinout1, int spinout2, int photonpol, 
				   int shellindex1, int shellindex2, int two_m1, int two_m2){

  static FourVector<complex<double> > polVectorPlus(0.,
						    -1./sqrt(2.),
						    complex<double>(0.,-1./sqrt(2.)),
						 0.);
  static FourVector<complex<double> > polVectorMin(0.,
						   1./sqrt(2.),
						   complex<double>(0.,-1./sqrt(2.)),
						   0.);
  static FourVector<complex<double> > polVector0(1.,0.,0.,0.);
  static FourVector<complex<double> > polVectorX(0.,1.,0.,0.);
  static FourVector<complex<double> > polVectorY(0.,0.,1.,0.);
  static FourVector<complex<double> > polVectorZ(0.,1.,0.,0.);
  pf1 = tk.GetN3();
  pf2 = tk.GetK3();
  qvec3 = tk.GetG3();
  
  J1=new NucleonEMOperator(tk.GetQsquared(),particletype1,0);
  J2=new NucleonEMOperator(tk.GetQsquared(),particletype2,0);
  FourVector<double> q(tk.GetWlab(),0.,0.,tk.GetKlab());
  GammaStructure Jcontr1, Jcontr2;
  if(photonpol==0){
    Jcontr1 = J1->getCC2(q, 0.,0,*pnucl)*polVector0;
    Jcontr2 = J2->getCC2(q, 0.,0,*pnucl)*polVector0;
  }
  else if(photonpol==1){
    Jcontr1 = J1->getCC2(q, 0.,0,*pnucl)*polVectorPlus;
    Jcontr2 = J2->getCC2(q, 0.,0,*pnucl)*polVectorPlus;
  }
  else if(photonpol==-1){
    Jcontr1 = J1->getCC2(q, 0.,0,*pnucl)*polVectorMin;
    Jcontr2 = J2->getCC2(q, 0.,0,*pnucl)*polVectorMin;
  }
  else{ cerr << "invalid photon pol" << endl;  exit(1); }
  FastParticle nucl1(particletype1, 0, pf1,tk.GetQsquared()/1.e06,0.,1.,1.,homedir);
  FastParticle nucl2(particletype2, 0, pf2,tk.GetQsquared()/1.e06,0.,1.,1.,homedir);
  if(!pw){
    gridf1 = new GlauberGridThick(60,18,5,pnucl,prec,integration_method(integrator),homedir);
    gridf2 = new GlauberGridThick(60,18,5,pnucl,prec,integration_method(integrator),homedir);
    gridf1->addParticle(nucl1);
    gridf2->addParticle(nucl2);
    //gridf1->printParticles();
    //gridf2->printParticles();
    gridf1->fillGrids();
    gridf2->fillGrids();
    gridf1->clearKnockout();
    gridf2->clearKnockout();
    gridf1->addKnockout(shellindex1,two_m1);
    gridf1->addKnockout(shellindex2,two_m2);
    gridf2->addKnockout(shellindex1,two_m1);
    gridf2->addKnockout(shellindex2,two_m2);
  }
  HelicitySpinor out1(tk.GetN4(),tk.GetMn(),(spinout1==-1?Spinor::Polarization::kDown:Spinor::Polarization::kUp));
  HelicitySpinor out2(tk.GetK4(),tk.GetMk(),(spinout2==-1?Spinor::Polarization::kDown:Spinor::Polarization::kUp));
  HelicitySpinor interm1down(tk.GetN4()-tk.GetG4(),tk.GetMn(),Spinor::Polarization::kUp);
  HelicitySpinor interm1up(tk.GetN4()-tk.GetG4(),tk.GetMn(),Spinor::Polarization::kDown);
  HelicitySpinor interm2down(tk.GetK4()-tk.GetG4(),tk.GetMk(),Spinor::Polarization::kUp);
  HelicitySpinor interm2up(tk.GetK4()-tk.GetG4(),tk.GetMk(),Spinor::Polarization::kDown);
  
  J1contrdown = Spinor::Bar(out1)*Jcontr1*interm1down;
  J1contrup = Spinor::Bar(out1)*Jcontr1*interm1up;
  J2contrdown = Spinor::Bar(out2)*Jcontr2*interm2down;
  J2contrup = Spinor::Bar(out2)*Jcontr2*interm2up;
  
  pair1up = new Pair("paris",homedir, pnucl->getA(),
		      pnucl->getN_array()[shellindex1],
		      pnucl->getL_array()[shellindex1],
		      1, pnucl->getJ_array()[shellindex1],
		      two_m1, -2*particletype1+1,
		      pnucl->getN_array()[shellindex2],
		      pnucl->getL_array()[shellindex2],
		      spinout2, pnucl->getJ_array()[shellindex2],
		      two_m2, -2*particletype2+1);
  pair1down = new Pair("paris",homedir, pnucl->getA(),
		      pnucl->getN_array()[shellindex1],
		      pnucl->getL_array()[shellindex1],
		      -1, pnucl->getJ_array()[shellindex1],
		      two_m1, -2*particletype1+1,
		      pnucl->getN_array()[shellindex2],
		      pnucl->getL_array()[shellindex2],
		      spinout2, pnucl->getJ_array()[shellindex2],
		      two_m2, -2*particletype2+1);
  pair2up = new Pair("paris",homedir, pnucl->getA(),
		      pnucl->getN_array()[shellindex1],
		      pnucl->getL_array()[shellindex1],
		      spinout1, pnucl->getJ_array()[shellindex1],
		      two_m1, -2*particletype1+1,
		      pnucl->getN_array()[shellindex2],
		      pnucl->getL_array()[shellindex2],
		      1, pnucl->getJ_array()[shellindex2],
		      two_m2, -2*particletype2+1);
  pair2down = new Pair("paris",homedir, pnucl->getA(),
		      pnucl->getN_array()[shellindex1],
		      pnucl->getL_array()[shellindex1],
		      spinout1, pnucl->getJ_array()[shellindex1],
		      two_m1, -2*particletype1+1,
		      pnucl->getN_array()[shellindex2],
		      pnucl->getL_array()[shellindex2],
		      -1, pnucl->getJ_array()[shellindex2],
		      two_m2, -2*particletype2+1);
		    
  //gsl MC setup			
  double res_real, res_imag, err_real, err_imag;

  double xl[6] = { 0., -1., 0., 0., -1., 0.};
  double xu[6] = { pnucl->getRange(), 1., 2.*PI, pnucl->getRange(), 1., 2.*PI };

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function real = { &matrixel_MC_real, 6, this};
  gsl_monte_function imag = { &matrixel_MC_imag, 6, this};

  size_t calls = 1.E6;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
  gsl_monte_vegas_init (s);
  s->stage=0;
  //s->verbose=1;
  gsl_monte_vegas_integrate (&real, xl, xu, 6, 1.E6, r, s,
			      &res_real, &err_real);
  display_results ("vegas warm-up real", res_real, err_real);

  cout << "converging..." << endl;

  s->stage=1;
  do
    {
      gsl_monte_vegas_integrate (&real, xl, xu, 6, calls, r, s,
				  &res_real, &err_real);
      printf ("result = % .6f sigma = % .6f "
	      "chisq/dof = %.1f\n", res_real, err_real, s->chisq);
      cout << flush;
      s->stage=3;
    }
//   while (fabs (s->chisq - 1.0) > 0.5);
  while (abs(res_real)>1.E-09? (abs(err_real/res_real)>1.E-02):0);

  
  
  gsl_monte_vegas_init (s);
  s->stage=0;
  //s->verbose=1;
  gsl_monte_vegas_integrate (&imag, xl, xu, 6, 1.E6, r, s,
			      &res_imag, &err_imag);
  display_results ("vegas warm-up imag", res_imag, err_imag);

  printf ("converging...\n");

  s->stage=1;
  do
    {
      gsl_monte_vegas_integrate (&imag, xl, xu, 6, calls, r, s,
				  &res_imag, &err_imag);
      printf ("result = % .6f sigma = % .6f "
	      "chisq/dof = %.1f\n", res_imag, err_imag, s->chisq);
      s->stage=3;
    }
//   while (fabs (s->chisq - 1.0) > 0.5);
  while (abs(res_imag)>1.E-09? (abs(err_imag/res_imag)>1.E-02):0);
  
  display_results ("vegas final imag", res_imag, err_imag);
  cout << abs(err_imag/res_imag) << endl << endl;
  gsl_monte_vegas_free (s);

  
 
  if(!pw){ delete gridf1; delete gridf2;}
  delete J1; delete J2;
  delete pair1up; delete pair1down; delete pair2up; delete pair2down;
//   complex<double> result;
//   if(CT) return results[1];
//   else return results[0];
  return complex<double>(res_real,res_imag);
}

complex<double> DoubleNModel::MC_helper(const double r1, const double costheta1, const double phi1, 
					const double r2, const double costheta2, const double phi2) const{
  
  double sintheta1 = sqrt(1.-costheta1*costheta1);
  double sintheta2 = sqrt(1.-costheta2*costheta2);
  
  nucqcd::Vector3 r1vec = nucqcd::Vector3(sintheta1*cos(phi1),sintheta1*sin(phi1),costheta1)*r1;
  nucqcd::Vector3 r2vec = nucqcd::Vector3(sintheta2*cos(phi2),sintheta2*sin(phi2),costheta2)*r2; // z was phi2 until 2026-11
  
  complex<double> glf1r1, glf1r2, glf2r1, glf2r2;
  if(getPW()){ glf1r1=glf1r2=glf2r1=glf2r2=1.; } // no FSI grids in plane-wave mode (this crashed until 2026-11)
  else if(getSRC()){
    if(getCT()){
      glf1r1=getGridf1()->getFsiSrcCtGridFull_interp3(r1,costheta1,phi1);
      glf1r2=getGridf1()->getFsiSrcCtGridFull_interp3(r2,costheta2,phi2);
      glf2r1=getGridf2()->getFsiSrcCtGridFull_interp3(r1,costheta1,phi1);
      glf2r2=getGridf2()->getFsiSrcCtGridFull_interp3(r2,costheta2,phi2);      
    }
    else{
      glf1r1=getGridf1()->getFsiSrcGridFull_interp3(r1,costheta1,phi1);
      glf1r2=getGridf1()->getFsiSrcGridFull_interp3(r2,costheta2,phi2);
      glf2r1=getGridf2()->getFsiSrcGridFull_interp3(r1,costheta1,phi1);
      glf2r2=getGridf2()->getFsiSrcGridFull_interp3(r2,costheta2,phi2);      
    }
  }
  else{
    if(getCT()){
      glf1r1=getGridf1()->getFsiCtGridFull_interp3(r1,costheta1,phi1);
      glf1r2=getGridf1()->getFsiCtGridFull_interp3(r2,costheta2,phi2);
      glf2r1=getGridf2()->getFsiCtGridFull_interp3(r1,costheta1,phi1);
      glf2r2=getGridf2()->getFsiCtGridFull_interp3(r2,costheta2,phi2);      
    }
    else{
      glf1r1=getGridf1()->getFsiGridFull_interp3(r1,costheta1,phi1);
      glf1r2=getGridf1()->getFsiGridFull_interp3(r2,costheta2,phi2);
      glf2r1=getGridf2()->getFsiGridFull_interp3(r1,costheta1,phi1);
      glf2r2=getGridf2()->getFsiGridFull_interp3(r2,costheta2,phi2);      
    }    
  }
  
  
  return ((getJ1contrup()*getPair1up()->getwf(getCorr(),r1,costheta1,phi1,r2,costheta2,phi2)+
		getJ1contrdown()*getPair1down()->getwf(getCorr(),r1,costheta1,phi1,r2,costheta2,phi2))
	       * exp(-nucqcd::i_unit*((getPf1()-getQvec3())*r1vec+getPf2()*r2vec)*nucqcd::inv_hbarc)
		* glf1r1*glf2r2
	- (getJ2contrup()*getPair2up()->getwf(getCorr(),r2,costheta2,phi2,r1,costheta1,phi1)+
		getJ2contrdown()*getPair2down()->getwf(getCorr(),r2,costheta2,phi2,r1,costheta1,phi1))
	       * exp(-nucqcd::i_unit*((getPf2()-getQvec3())*r1vec+getPf1()*r2vec)*nucqcd::inv_hbarc)
		* glf1r2*glf2r1)*r1*r1*r2*r2;
  
  
  
}
  

//x[0] r1, x[1] costh1, x[2] phi1, x[3]-x[5] similar for vec{r_2}
//p is pointer to DoubleNmodel object.
//dim has to be 6!
double matrixel_MC_real (double* x, size_t dim, void * p){
  
  if(dim!=6){
    cerr << "Dim of function in MC has to be 6!" << endl;
    exit(1);
  }
  DoubleNModel *model = (DoubleNModel *)p;
  return real(model->MC_helper(x[0],x[1],x[2],x[3],x[4],x[5]));
}

double matrixel_MC_imag (double* x, size_t dim, void * p){
  
  if(dim!=6){
    cerr << "Dim of function in MC has to be 6!" << endl;
    exit(1);
  }
  DoubleNModel *model = (DoubleNModel *)p;
  return imag(model->MC_helper(x[0],x[1],x[2],x[3],x[4],x[5]));

}

void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  cout << flush;
  }

