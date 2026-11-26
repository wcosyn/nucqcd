#include <Tools/Core/IntegratorFlag.hpp>
#include <Knockout/RhoTCross.hpp>
#include <Tools/Core/AuxFunction.hpp>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/Utilfunctions.hpp>

namespace nucqcd {


using namespace std;

RhoTCross::RhoTCross(const int nucleus, const double p_max, const std::string dir, const bool no_cuts,
  const bool user_set, const double user_sigma, const double precision, const int integr, const int max_Eval,
  const double p_dil, const int particletype
):
homedir(dir),
pmax(p_max),
nocuts(no_cuts),
userset(user_set),
usersigma(user_sigma),
nucleusthick(nucleus,dir),
prec(precision),
integrator(nucqcd::checked_integrator_flag(integr)),
abserror(1.E-012),
maxEval(max_Eval),
pdil(p_dil),
particle_type(particletype){

//   pdistgrid = new DistMomDistrGrid*[nucleusthick.getTotalLevels()];
  for(int i=0;i<nucleusthick.getTotalLevels();i++){
    pfsigrid.push_back(std::make_unique<GlauberDecayGridThick>(60,20,5,&nucleusthick,prec,2,homedir));
//     pdistgrid[i] = new DistMomDistrGrid(i, pmax, 30,20,5,pfsigrid[i],1.E-03,2,2E04,0.,homedir);
  }
  nrofcross=pfsigrid[0]->getNumber_of_grids()+1;
}


//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoTCross::getCrosst(double *results, const double Ebeam, const double Q2, const double nu, const double t){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
//   double A = (t+Q2-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06)*0.5;
//   double MA=getNucleusthick().getMassA()*1.E-03;
//   double C0=nu+MA;
//   double pzrho=(MA*MA+qvec*qvec-C0*C0-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-2.*C0*A/nu)/(2.*qvec*(1.-C0/nu));
//   double Erho = (pzrho*qvec-A)/nu;
//   double pxrho = sqrt(Erho*Erho-pzrho*pzrho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06); //pxrho defines x-axis (so only positive sqrt is considered)
//   double prho = sqrt(Erho*Erho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06); //pxrho defines x-axis (so only positive sqrt is considered)
//   if(std::isnan(pxrho)) pxrho=0.;
//   double pAz=qvec-pzrho;
//   double z=Erho/nu;
//   if(z>1.||z<0.9){for(int i=0;i<nrofcross;i++) results[i]=0.;}
//   //else for(int i=0;i<nrofcross;i++){ results[i]=1./(2.*qvec*abs(C0-Erho+Erho*(1.-qvec*pzrho/prho/prho)));}
//   else { for(int i=0;i<nrofcross;i++) results[i]=1./(2.*MA);}
//   return;
  
  
  int res=90;
  unsigned count=0;
  if (integrator==1||integrator==2){
    numint::vector_d ret(nrofcross,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pmax*1.E-03,1.,2.*PI}};
    
    RhoTCross::Ftor_rho F;
    F.cross = this;
    F.Q2=Q2;
    F.nu=nu;
    F.qvec=qvec;
    F.t=t;
    F.Erho=0.;
    F.prho=0.;
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_rho::exec;
    mdf.param = &F;
    F.f=klaas_rho_t;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,1E03,maxEval,ret,count,0);
    for(int i=0;i<nrofcross;++i) results[i]=ret[i];
    if(ret[2]*1.E-04>abserror) {abserror=ret[2]*1.E-04; }
  }
  for(int i=0;i<nrofcross;i++) results[i]*= nucqcd::alpha_em*(Ebeam-nu)/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon))*pow(nucqcd::inv_hbarc*1.E03,3.)/nucleusthick.getA()/(2.*qvec);
//   cout << nu << " " << t << " " << results[0] << " " << res << " " << count << " " << abserror << " " << endl;
}




//input in GeV!!!!
//cross section is obtained in units GeV.  dsigma/dt(t=0) is unspecified in the formula (dimensions energy^-4).
void RhoTCross::getCrossz(double *results, const double Ebeam,  const double Q2, const double nu, const double z){

  double qvec = sqrt(Q2+nu*nu);
  double fScatterAngle = 1. - Q2/2./Ebeam/(Ebeam-nu);
  double epsilon = 1./(1.+2.*(1.+nu*nu/Q2)
		     *((1.-fScatterAngle)/(1.+fScatterAngle)));
  
  double Erho=nu*z;
  double prho=sqrt(Erho*Erho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06);
  
//   double MA=getNucleusthick().getMassA()*1.E-03;
//   double EA=-Erho+MA+nu;
//   double pA=sqrt(EA*EA-MA*MA);
//   if(std::isnan(pA)) pA=0.;
//   double pzrho = (prho*prho+qvec*qvec-pA*pA)/(2.*qvec);
//   double pxrho=sqrt(prho*prho-pzrho*pzrho);
//   if(std::isnan(pxrho)) pxrho=0.;
//   double pxA=-pxrho;
//   double pzA=qvec-pzrho;
//   double t = -Q2+nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-2.*Erho*nu+2.*qvec*pzrho;
//   if(t>-0.1||t<-0.4){for(int i=0;i<nrofcross;i++) results[i]=0.;}
//   else for(int i=0;i<nrofcross;i++){ results[i]=nu;}
//   return;
  
  int res=90;
  unsigned count=0;
  if (integrator==1||integrator==2){
    numint::vector_d ret(nrofcross,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{pmax*1.E-03,1.,2.*PI}};
    
    RhoTCross::Ftor_rho F;
    F.cross = this;
    F.Q2=Q2;
    F.nu=nu;
    F.qvec=qvec;
    F.t=0.;
    F.Erho=Erho;
    F.prho=prho;
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_rho::exec;
    mdf.param = &F;
    F.f=klaas_rho_z;
    if(integrator==1) res = numint::cube_romb(mdf,lower,upper,abserror,1.E-03,ret,count,0);
    else res = numint::cube_adaptive(mdf,lower,upper,abserror,1.E-03,1E03,maxEval,ret,count,0);
    for(int i=0;i<nrofcross;++i) results[i]=ret[i];
    if(ret[2]*1.E-04>abserror)  {abserror=ret[2]*1.E-04; }
  }
  
  for(int i=0;i<nrofcross;i++) results[i]*= nucqcd::alpha_em*(Ebeam-nu)/(2.*pow(2.*PI,2.)*Ebeam*Q2*(1.-epsilon))*pow(nucqcd::inv_hbarc*1.E03,3.)/nucleusthick.getA()*nu/qvec;
//   cout << nu << " " << z << " " << results[0] << " " << res << " " << count << " " << abserror << endl;
  
}




//all energies in MeV please
void RhoTCross::getMomdistr(double *results, double prho, double thetarho, double Q2, int shell, 
			    double pm, double pmcostheta, double pmphi){
  //if userset there's only one possible glaubergrid, so we don't have to check every time!
  if(!(userset&&pfsigrid[shell]->getFilledallgrid())){
    FastParticle rho(particle_type, 0, prho,0.,0.,Q2,145.,1.,1.,homedir);
    if(userset) rho.setScatter(usersigma,6.,-0.2,pdil);
    pfsigrid[shell]->clearParticles();
    pfsigrid[shell]->addParticle(rho);
    pfsigrid[shell]->updateGrids();
    pfsigrid[shell]->clearKnockout();
  }
  //update the grid if necessary
  nucqcd::Rotation rot;
  rot.Rotate(-thetarho,nucqcd::Vector3(0.,1.,0.));
  double key=double(shell)/getNucleusthick().getTotalLevels()+int(round(thetarho*10.));
  map<double,DistMomDistrGrid>::iterator it=distgridmap.find(key);
//   cout << key << " " << thetarho << " " << shell  << endl;
  if(it==distgridmap.end()){
    distgridmap[key]=DistMomDistrGrid(shell, pmax, 30,20,5,pfsigrid[shell].get(),1.E-03,2,5E04,0.,homedir);
//     distgridmap.insert(pair<double,DistMomDistrGrid>(key,DistMomDistrGrid(shell, pmax, 30,20,5,pfsigrid[shell].get(),1.E-03,2,2E04,0.,homedir)));
    distgridmap[key].updateGrids(pfsigrid[shell].get(),shell,rot);
//     distgridmap[key].printRho_grid(0);
//     cout << endl << endl;
    it=distgridmap.find(key);
  }
  else distgridmap[key].updateGrids(pfsigrid[shell].get(),shell,rot);

  results[nrofcross-1] = distgridmap[key].getRhopwGridFull_interp(pm);
  //fsi
  for(int i=0;i<nrofcross-1;i++) results[i] = distgridmap[key].getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
//    if(pm>200.&&results[0]>50.) {distgridmap[key].printRho_grid(0); cout  << endl << endl << endl;}
  
  
//   pdistgrid[shell]->updateGrids(pfsigrid[shell].get(),shell,rot);
//   //plane-wave
//   results[nrofcross-1] = pdistgrid[shell]->getRhopwGridFull_interp(pm);
//   //fsi
//   for(int i=0;i<nrofcross-1;i++) results[i]= pdistgrid[shell]->getRhoGridFull_interp3(i, pm, pmcostheta, pmphi);
}
  
  
double RhoTCross::getfrontfactor(double nu, double qvec, double Erho, double prho, double pzrho, double pxrho,
				 double s, double Q2, double mN, double t, bool torz, double C0, double Cz, double Cx){
  //elementary rho production cross section parametrized as exponential e^{beta*t}
  double elementary = (s*s-2.*s*(mN*mN-Q2)+pow(mN*mN+Q2,2.))/(mN*mN)*exp(t*6.);
  if(torz) return elementary;
  //z integration: the integrand has to carry the Jacobian K = |dt/dz|/(2 nu) between the two ways of eliminating the
  //on-shell delta function of the final nucleon, f(z,t) = (C-p_rho)^2-mN^2 = 0 with C = (C0,Cx,0,Cz) the four-momentum
  //of the struck nucleon plus photon, at fixed missing momentum:
  //  K = q |C0 px-Cz nu px/q-Cx nu (z-pz/q)| / |Cz px-Cx pz|.
  //For pm=0 this is C0-nu = massi, the off-shell energy M_A-M_(A-1)-E_exc of the struck nucleon (not the on-shell mN).
  //Derivation and checks: docs/rho_jacobian.md.  (Until 2026-11 the mass of the residual nucleus system, ~M_A, stood
  //here, which made the z integrals a factor ~A too large compared to the t integrals.)
  double z = Erho/nu;
  double den = std::abs(Cz*pxrho-Cx*pzrho);
  if(den<1.E-12) return 0.;
  double jac = qvec*std::abs(C0*pxrho-Cz*nu*pxrho/qvec-Cx*nu*(z-pzrho/qvec))/den;
  return jac*elementary;
}

void RhoTCross::klaas_rho_t(numint::vector_d & results, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
    double t, double dummy, double dummy2){
  results=numint::vector_d(cross.getNrofcross(),0.);
  if(pm==0.){
    return;
  }
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  double sintheta = sqrt(1.-costheta*costheta);
  
  
  //ALL IN GeV to avoid some numeric overflow almost zero shit!
  for(int i=0;i<cross.getNucleusthick().getTotalLevels();i++){
    double massi = cross.getNucleusthick().getMassA()-cross.getNucleusthick().getMassA_min_1(i)-cross.getNucleusthick().getExcitation()[i];
    massi*=1.E-03;
    double mN = i<cross.getNucleusthick().getPLevels()? nucqcd::mass_p:nucqcd::mass_n;
    mN*=1.E-03;
    //determine kinematics
    double A = (t+Q2-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06)*0.5;
    double C0 = nu+massi;
    double Cz = qvec+pm*costheta;
    double Cy = pm*sintheta*sinphi;
    double Cx = pm*sintheta*cosphi;
    double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam 
    double D = (s-mN*mN+nucqcd::mass_rho*nucqcd::mass_rho*1.E-06)/2.+C0*A/nu;
    double E = Cz-C0*qvec/nu;
    double a = E*E-Q2*Cx*Cx/(nu*nu);
    double b = 2.*(D*E+A*qvec*Cx*Cx/(nu*nu));
    double c = D*D-A*A*Cx*Cx/(nu*nu)+Cx*Cx*nucqcd::mass_rho*nucqcd::mass_rho*1.E-06;
    double discr = b*b-4.*a*c;
    if(discr>-1.E-09){
      //zero discriminant
      if(abs(discr)<1.E-09){
	double pzrho=-b/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	//check cuts
	if(cross.getNocuts()||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(std::isnan(pxrho)) pxrho=0.; //underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);	  
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1,C0,Cz,Cx);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double Erho = (pzrho*qvec-A)/nu;
	//check cuts
	if(cross.getNocuts()||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(std::isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1,C0,Cz,Cx);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	Erho = (pzrho*qvec-A)/nu;
	//check cuts
	if(cross.getNocuts()||Erho/nu>0.9){
	  double pxrho = sqrt(Erho*Erho-pzrho*pzrho-nucqcd::mass_rho*nucqcd::mass_rho*1.E-06);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(std::isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(D+E*pzrho)==SIGN(-Cx*pxrho))||pxrho==0.){ 
	    double prho = sqrt(pzrho*pzrho+pxrho*pxrho);
	    double intresults[cross.getNrofcross()];	  
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,1,C0,Cz,Cx);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }    
    }
  }
  for(int i=0;i<cross.getNrofcross();i++){
    results[i]*=pm*pm;
  }
}

void RhoTCross::klaas_rho_z(numint::vector_d & results, double pm, double costheta, double phi, RhoTCross & cross, double Q2, double nu, double qvec,
    double dummy, double Erho, double prho){
  results=numint::vector_d(cross.getNrofcross(),0.);
  if(pm==0.){
    return;
  }
  double sintheta = sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  
  for(int i=0;i<cross.getNucleusthick().getTotalLevels();i++){
    //offshell initial energy
    double massi = cross.getNucleusthick().getMassA()-cross.getNucleusthick().getMassA_min_1(i)-cross.getNucleusthick().getExcitation()[i];
    massi*=1.E-03;
    double mN = i<cross.getNucleusthick().getPLevels()? nucqcd::mass_p:nucqcd::mass_n;
    mN*=1.E-03;
    double C0 = nu+massi;
    double Cz = qvec+pm*costheta;
    double Cy = pm*sintheta*sinphi;
    double Cx = pm*sintheta*cosphi;
    double s = C0*C0-Cz*Cz-Cx*Cx-Cy*Cy; //mandelstam 
    double En = C0-Erho;
    double pn = sqrt(En*En-mN*mN);
    
    double A=(Cx*Cx+Cy*Cy+Cz*Cz-pn*pn+prho*prho)*0.5;
    double a = Cz*Cz+Cx*Cx;
    double b = -2.*A*Cz;
    double c = A*A-Cx*Cx*prho*prho;
    double discr = b*b-4.*a*c;
    if(discr>-1.E-09){
      if(abs(discr)<1.E-09) {
	double pzrho = -b/(2.*a);
	double t =  -Q2 + nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(cross.getNocuts()||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)           
	  if(std::isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0,C0,Cz,Cx);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }
      else{
	discr = sqrt(discr);
	double pzrho = (-b+discr)/(2.*a);
	double t =  -Q2 + nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(cross.getNocuts()||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(std::isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 

	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0,C0,Cz,Cx);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
	pzrho = (-b-discr)/(2.*a);
	t =  -Q2 + nucqcd::mass_rho*nucqcd::mass_rho*1.E-06-2*Erho*nu+2.*qvec*pzrho;
	//check cuts
	if(cross.getNocuts()||(t<-0.1&&t>-0.4)){
	  double pxrho = sqrt(prho*prho-pzrho*pzrho);//pxrho defines x-axis (so only positive sqrt is considered)     
	  if(std::isnan(pxrho)) pxrho=0.;//underflow
	  if((SIGN(A-Cz*pzrho)==SIGN(Cx*pxrho))||pxrho==0.){ 
	  
	    double intresults[cross.getNrofcross()];
	    cross.getMomdistr(intresults,prho*1.E03,atan2(pxrho,pzrho),Q2,i,pm*1.E03,costheta,phi);
	    double front=cross.getfrontfactor(nu,qvec,Erho,prho,pzrho,pxrho,s,Q2,mN,t,0,C0,Cz,Cx);
	    //correct for not completely full shells
	    if(i==cross.getNucleusthick().getPLevels()-1) front*=double((cross.getNucleusthick().getFinalMProton()+1)*2-cross.getNucleusthick().getOnlyOneProton())/(cross.getNucleusthick().getJ_array()[i]+1);
	    if(i==cross.getNucleusthick().getTotalLevels()-1) front*=double((cross.getNucleusthick().getFinalMNeutron()+1)*2-cross.getNucleusthick().getOnlyOneNeutron())/(cross.getNucleusthick().getJ_array()[i]+1);
	    for(int dd=0;dd<cross.getNrofcross();dd++) results[dd]+=front*intresults[dd];
	  }
	}
      }
    }
  }
  
  for(int i=0;i<cross.getNrofcross();i++){
    results[i]*=pm*pm;
  }
}

}  // namespace nucqcd

  