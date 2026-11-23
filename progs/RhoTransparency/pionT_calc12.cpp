#include <iostream>
#include <cstdlib>

using namespace std;

#include <Knockout/PionTCross.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/RunLog.hpp>
#include <vector>
#include <string>



int main(int argc, char *argv[])
{

    const std::vector<std::string> arg_names={"nucleus name", "beam energy [GeV]", "scattered energy [GeV]", "scattered angle [deg]"};
    if(argc!=5){
        cout << "Usage: " << argv[0] << " <nucleus name> <beam energy [GeV]> <scattered energy [GeV]> <scattered angle [deg]>" << endl;
        return 1;
    }
    nucqcd::RunLog runlog(argc, argv, arg_names);  // run bookkeeping as "# " lines (Tools/Core/RunLog.hpp)
    std::cout << "# compiled from file: " << __FILE__ << std::endl;



    string nucleus_name=argv[1];  //nucleus name: options are He, C, O, Fe, Pb, Al, Cu, Au
    if(nucleus_name!="He" && nucleus_name!="C" && nucleus_name!="O" && nucleus_name!="Fe" && nucleus_name!="Pb" && nucleus_name!="Al" && nucleus_name!="Cu" && nucleus_name!="Au"){
        cout << "Nucleus name not recognized. Options are He, C, O, Fe, Pb, Al, Cu, Au." << endl;
        return 1;
    }
    double Ebeam  = atof(argv[2])*1.E03;  //Beam energy in GeV, converted to MeV
    double E_out = atof(argv[3])*1.E03;  // Scattered electron energy in GeV, converted to MeV
    double theta_e = atof(argv[4])*nucqcd::deg_to_rad; // Scattered electron angle in degrees input converted to radians
    
    double Q2 = 2.*Ebeam*E_out*(1.-cos(theta_e)); //Q2 in MeV^2

    PionTCross pionT_object = PionTCross(MeanFieldNucleus::TypeNames.at(nucleus_name),
                                  400.,nucqcd::share_dir,0,0.,1.E-05,2,20000);
    double *results = new double[pionT_object.getNrofcross()];
    pionT_object.getCross(results, Ebeam,E_out,theta_e);

    cout << Q2*1.E-06 << " ";
    for(int i=0;i<pionT_object.getNrofcross();i++){
        cout << results[i] << " ";
    }
    for(int i=0;i<pionT_object.getNrofcross()-1;i++){
        cout << results[i]/results[pionT_object.getNrofcross()-1] << " ";
    }
    cout << endl;

    delete [] results;
    return 0;
}