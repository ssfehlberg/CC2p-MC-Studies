#include "constants.h"
using namespace Constants;
class helper_funcs{
 
 public:
  virtual void Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual bool In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z);
  virtual double real_sqrt(double x);
  virtual double Which_Energy(bool total_E, TVector3 mom, double mass);

  bool fv; //is the overlay event in the FV?

};//end of class definition


//Used in the Overlay to help with the MC breakdown definitions. Does same thing as above but returns fv variable
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void helper_funcs::Overlay_In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
  float_t xmin = 0.0 + x_low_edge;
  float_t xmax = 256.35 - x_up_edge;
  float_t ymin = -116.5 + y_low_edge;
  float_t ymax = 116.5 - y_up_edge;
  float_t zmin = 0.0 + z_low_edge;
  float_t zmax = 1036.8 - z_up_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    fv = false;
  } else{
    fv = true;
  } 
}

//Determines if the input x,y,z is wihtin the FV. The distance from any edge can be defined by VAR_low(high)_edge
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool helper_funcs::In_FV(float x_low_edge, float x_up_edge, float y_low_edge, float y_up_edge, float z_low_edge, float z_up_edge, float x, float y, float z){
  float_t xmin = 0.0 + x_low_edge;
  float_t xmax = 256.35 - x_up_edge;
  float_t ymin = -116.5 + y_low_edge;
  float_t ymax = 116.5 - y_up_edge;
  float_t zmin = 0.0 + z_low_edge;
  float_t zmax = 1036.8 - z_up_edge;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    return false;
  } else{
    return true;
  } 
}



//Defines Real Square Root in Case number is negative
//////////////////////////////////////////////////////
double helper_funcs::real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

//Calculates either the kinetic or total energy of a particle. Total_e = true -> total energy, Total_E = false-> KE. Needs momentum and mass of particle                                                                                          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double helper_funcs::Which_Energy(bool total_E, TVector3 mom, double mass){
  double Energy;
  if(total_E == true){ //Total Energy
    Energy = std::sqrt(mom.Mag2() + std::pow(mass,2));
  } else if (total_E == false){
    Energy = std::sqrt(mom.Mag2() + std::pow(mass,2)) - mass; //Kinetic Energy
  }
  return Energy;
}

