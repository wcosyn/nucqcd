/*
 * FsiCorrelator_Solver.cpp
 * FsiCorrelator: construction of the correlation grid
 * (the matrix of the equations, the correlation function and the Newton solver).
 * Part of the implementation of FsiCorrelator, split from FsiCorrelator.cpp in 2026-11 (the functions were moved verbatim).
 */

#include <fstream>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <Glauber/FsiCorrelator.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/GridCache.hpp>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {

using namespace std;

#define TOTALGRIDPNTS (getCorr_rgrid())*(getCorr_cthgrid())

void FsiCorrelator::constructCorrgrid(const double *density, std::vector<std::vector<double>> &corrgrid){

  cout << "Constructing correlation grid..." << endl;
  functionmatrix.assign(TOTALGRIDPNTS,std::vector<double>(TOTALGRIDPNTS));//TOTALGRIDPNTS= getCorr_rgrid()*getCorr_thgrid()
  cout << "Construct Matrix of the functions..." << endl;
  constructMatrixSet(density);
  x.assign(TOTALGRIDPNTS,1.);
  int check=0;
  newt(&check);
  cout << ReturnFunction() << endl;
  functionmatrix.clear();
  if(check>0) {
    x.clear();
    throw std::runtime_error(nucqcd::ErrorMessage() << "no solution found for integration equation" );
  }
  else{
    for(int i=0;i<TOTALGRIDPNTS;i++){
      int j = i%(corr_cthgrid);
      int ii = i/(corr_cthgrid);
      //cout << i << " " << ii << " " << j << " " << x[i] << endl;
      corrgrid[ii][j]=x[i];
    }
    cout << endl << endl;
    x.clear();
    return;
  }  
}



//////////////////////
//construct the matrix needed to solve the non-linear set of equations
////////////////////////

void FsiCorrelator::constructMatrixSet(const double *density){

  double VOLEL=pnucleus->getRange()/corr_rgrid*2./corr_cthgrid*PI/corr_phigrid; //phasespace
  //create grid, phi dependence is ignored as it is "very small"
  for(int rindex1=0; rindex1<corr_rgrid; rindex1++){
    double r1 = pnucleus->getRange()*rindex1/corr_rgrid;
    for(int thetaindex1=0; thetaindex1<corr_cthgrid; thetaindex1++){
      //double theta1 = PI*thetaindex1/corr_thgrid;
      double costheta1 = 1.-thetaindex1/invcthstep;
      double sintheta1 = sqrt(1.-costheta1*costheta1);
/*      double costheta1 = cos(theta1);
      double sintheta1 = sin(theta1);*/
//       double costheta1,sintheta1;
//       sincos(theta1,&sintheta1,&costheta1);
      int index1 = (rindex1)*(corr_cthgrid)+thetaindex1;
      for(int rindex2=0; rindex2<corr_rgrid; rindex2++){
	double r2 = pnucleus->getRange()*rindex2/corr_rgrid;
	for(int thetaindex2=0; thetaindex2<corr_cthgrid; thetaindex2++){
	  double costheta2 = 1.-thetaindex2/invcthstep;
	  double sintheta2 = sqrt(1.-costheta2*costheta2);
// 	  double theta2 = PI*thetaindex2/corr_thgrid;
// 	  double costheta2,sintheta2;
// 	  sincos(theta2,&sintheta2,&costheta2);
	  int index2 = (rindex2)*(corr_cthgrid)+thetaindex2;	  
	  functionmatrix[index2][index1] = pnucleus->getDensity(r1,density)*VOLEL;
	  double gsum = 0.;
	  for(int phiindex1=0; phiindex1<2*corr_phigrid; phiindex1++){
	    double phi1 = PI*phiindex1/corr_phigrid;
	    double cosphi1, sinphi1;
	    sincos(phi1,&sinphi1,&cosphi1);
	    for(int phiindex2=0; phiindex2<2*corr_phigrid; phiindex2++){
	      double phi2 = PI*phiindex2/corr_phigrid;
	      double cosphi2, sinphi2;
	      sincos(phi2,&sinphi2,&cosphi2);
	      gsum +=correlation(normr(r1, costheta1, sintheta1, cosphi1, sinphi1, r2, costheta2, sintheta2, cosphi2, sinphi2));
		//begin and endpoint trapezium integration coefficients
		/*if(phiindex1==0) fillmatrix[index1][index2]/=2.;		
		if(phiindex1==2.*PHIGRID) fillmatrix[index1][index2]/=2.;
		if(thetaindex1==0) fillmatrix[index1][index2]/=2.;
		if(thetaindex1==THETAGRID) fillmatrix[index1][index2]/=2.;
		if(rindex1==0) fillmatrix[index1][index2]/=2.;
		if(rindex1==RGRID)  fillmatrix[index1][index2]/=2.;*/
	      }
	  }
	  functionmatrix[index2][index1] *=gsum;
// 	  if(thetaindex1==0) functionmatrix[index1][index2]/=2.;
// 	  if(thetaindex1==corr_cthgrid+1) functionmatrix[index1][index2]/=2.;
// 	  if(rindex1==0) functionmatrix[index1][index2]/=2.;
// 	  if(rindex1==corr_rgrid+1)  functionmatrix[index1][index2]/=2.;
	}
      }
    }
  }
}



/////////////////
//G-D correlation function
///////////////////////
  
double FsiCorrelator::correlation(double r) const{

  double a[11]={0.99994,0.99600E-01,-2.9385,28.513,-132.53,238.33,-151.54,-69.845,159.94,-88.043,16.930};
  if(r>1.24) return 1+0.1204071*exp(-2.5002*(r-1.24)*(r-1.24));
  else{
    double result=0.;
    for(int i=0;i<11;i++) result += a[i]*pow(r,i);
    return 1-result;
    }
  //return 1-exp(-r*r*1.18147448);
}

////
//some more num recipes shit i copied
///////


void FsiCorrelator::newt(int *check){

  int MAXITS=200;
  double TOLF = 1.0e-7, TOLMIN = 1.0e-06, TOLX = 1.0e-07, STPMX = 100.0;

  std::vector<double> functionvalues(TOTALGRIDPNTS);
  double function = ReturnFunction(functionvalues);
  double test=0.;
  //test if initial guess is zero point
  for(int i=0;i<TOTALGRIDPNTS;i++) if(abs(functionvalues[i])>test) test=abs(functionvalues[i]);
  if(test<TOLF) {
    cout << "initial point is zero " << test << endl;
    return;
  }
 
  double sum = 0.;
  for(int i=0;i<TOTALGRIDPNTS;i++) sum+=x[i]*x[i];
  double stpmax = STPMX*(sqrt(sum)>float(TOTALGRIDPNTS)? sqrt(sum):float(TOTALGRIDPNTS));
  for(int its=0;its<MAXITS;its++){
    cout << "Entering iteration " << its << endl;
    std::vector<std::vector<double> > Jacobian(TOTALGRIDPNTS,std::vector<double>(TOTALGRIDPNTS));
    std::vector<double> grad(TOTALGRIDPNTS);
    std::vector<double> xold(TOTALGRIDPNTS);
    double functionold;
    std::vector<double> step(TOTALGRIDPNTS);
    cout << "Constructing Jacobian...." << endl;
    ConstructJacobian(Jacobian);
    
    for(int i=0; i<TOTALGRIDPNTS; i++){
      sum=0.;
      for(int j=0;j<TOTALGRIDPNTS;j++) sum+=functionvalues[j]*Jacobian[j][i];
      grad[i] = sum;
    }
    for (int i=0;i<TOTALGRIDPNTS;i++) {xold[i] = x[i]; step[i] = -functionvalues[i];}
    functionold = function;
    cout << "Solving for step... " << TOTALGRIDPNTS << endl;
    LUsolve(Jacobian,step);
    // for(int i=0;i<TOTALGRIDPNTS;i++) cout << step[i] << " ";
    cout << "Doing linesearch..." << endl;
    lnsrch(xold,functionold,grad,step,&function,stpmax,check,functionvalues);
    test=0.;
    for(int i=0;i<TOTALGRIDPNTS;i++) if(abs(functionvalues[i])>test) test = abs(functionvalues[i]);
    if(test<TOLF) {
      cout << "Convergence reached on function values " << endl;
      *check=0;
      return;
    }
    
    if(*check){  //check for grad f equal to zero
      test=0.;
      double den=(function>0.5*TOTALGRIDPNTS? function : 0.5*TOTALGRIDPNTS);
      for(int i=0;i<TOTALGRIDPNTS;i++){
	double temp=abs(grad[i]) * (abs(x[i])>1.? abs(x[i]) : 1.)/den;
	if(temp>test) test=temp;
      }
      *check=(test < TOLMIN ? 1:0);
      cout << "grad f was equal to zero, spurious thingy" << endl;
      return;
    }
    test=0.;
    for(int i=0;i<TOTALGRIDPNTS;i++){
      double temp=(abs(x[i]-xold[i]))/(abs(x[i])>1.? abs(x[i]) : 1.);
      if(temp>test) test=temp;
    }
    if(test < TOLX) {
      cout << "Convergence reached on dx" << endl;
      return;
    }

  }
  
  cout << "Max iterations reached" << endl;
  *check = 2;
  return;

}

/////////////////
//Jacobian is contructed, what did i need it for again? :)
/////////////////

void FsiCorrelator::ConstructJacobian(std::vector<std::vector<double> > &Jacobian){

  for(int i=0;i<TOTALGRIDPNTS;i++){
    for(int j=0;j<TOTALGRIDPNTS;j++){
      Jacobian[i][j]=functionmatrix[i][j]*x[i];
      if(i==j) {
	for(int k=0;k<TOTALGRIDPNTS;k++) Jacobian[i][i]+=functionmatrix[i][k]*x[k];
      }
    }
  }
}


///////////////
//i think this is another matrix needed to solve the non-lin set of equations :)
//aha, i think i remember, this is the function that needs minimizing!
/////////////////

double FsiCorrelator::ReturnFunction(std::vector<double> &functionvalues){
  for(int i=0;i<TOTALGRIDPNTS;i++){
    functionvalues[i]=-1.*2*corr_phigrid;
    for(int j=0;j<TOTALGRIDPNTS;j++){
      functionvalues[i]+=functionmatrix[i][j]*x[i]*x[j];
    }
  }
  double f=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) f+=functionvalues[i]*functionvalues[i];
  f/=2.;
  return f;
}

double FsiCorrelator::ReturnFunction(){
  std::vector<double> functionvalues(TOTALGRIDPNTS);
  for(int i=0;i<TOTALGRIDPNTS;i++){
    functionvalues[i]=-1.*2*corr_phigrid;
    for(int j=0;j<TOTALGRIDPNTS;j++){
      functionvalues[i]+=functionmatrix[i][j]*x[i]*x[j];
    }
  }
  double f=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) f+=functionvalues[i]*functionvalues[i];
  f/=2.;
  return f;
}

///////////////
//taken from numerical recipes, needed in the procedure to solve non-lin set, by minimalizing a testfunction
///////////////////////

void FsiCorrelator::lnsrch(const std::vector<double> &xold, double functionold, const std::vector<double> &grad, std::vector<double> &step, double *functionnew,
	    double stpmax, int *check, std::vector<double> &functionvalues){

  double TOLX = 1.0e-07;
  double ALF = 1.0e-04; 

  *check=0;
  float sum=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) sum += step[i]*step[i];
  sum=sqrt(sum);
  if(sum > stpmax) for(int i=0;i<TOTALGRIDPNTS;i++) step[i] *= stpmax/sum;  //scaling
  
  double slope=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) slope += grad[i]*step[i];

  double test=0.;
  for(int i=0;i<TOTALGRIDPNTS;i++){
    double temp=abs(step[i])/(fabs(xold[i])>1.0? xold[i]:1.0);
    if (temp>test) test=temp;
  }
  double alamin=TOLX/test;  //compute minimum value of lambda
  double alam=1.0;  //start with full step

  for(;;){

    double alam2=0., functionnew2=0., functionold2=0., tmplam=0.;

    for(int i=0;i<TOTALGRIDPNTS;i++) x[i]=xold[i]+alam*step[i];
    *functionnew = ReturnFunction(functionvalues);
    if(alam <alamin) {
      for(int i=0;i<TOTALGRIDPNTS;i++) x[i] = xold[i];
      *check=1;
      return;
    }
    else if (*functionnew <= functionold + ALF*alam*slope) return;
    else{

      if(alam == 1.) tmplam = -slope/(2.0*(*functionnew-functionold-slope));
      else{
	double rhs1 = *functionnew-functionold-alam*slope;
	double rhs2 = functionnew2-functionold2-alam2*slope;
	double a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	double b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a==0.) tmplam = -slope/(2.0*b);
	else{
	  double disc = b*b-3.0*a*slope;
	  if(disc<0.) cout << "roundoff problem in lnsrch" << endl;
	  else tmplam = (-b+sqrt(disc))/(3.0*a);
	}
	if(tmplam>0.5*alam) tmplam = 0.5*alam;
      }
    }
    alam2=alam;
    functionnew2=*functionnew;
    functionold2=functionold;
    alam=(tmplam>0.1*alam? tmplam : 0.1*alam);
  }
}

}  // namespace nucqcd



