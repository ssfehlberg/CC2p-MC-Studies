#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants{

  //debug statements
  bool _debug = false; 

  //input and output directory
  const char* input = "/uboone/data/users/sfehlber/CC2p/MC_Studies/GENIE/samples";
  const char* output = "/uboone/data/users/sfehlber/CC2p/MC_Studies";
    
  //Are we adding protons together in the STVs?
  bool add_protons = true;

  //Are we using the xsecc binning
  bool use_xsec_binning = true;

  //FOR GENIE ONLY
  //do we want with FSI or without FSI?
  bool fsi = true;

  //Momentum Cut Values (only used by the overlay)
  double MUON_MOM_CUT_LOW = 0.1;
  double MUON_MOM_CUT_HIGH = 1.2;
  double PROTON_MOM_CUT_LOW = 0.3;
  double PROTON_MOM_CUT_HIGH = 1.0;
  double CHARGED_PI_MOM_CUT = 0.065;

  //Argon 40 properties
  double MASS_TARGET = 37.215526; //GeV
  // This binding energy value is used in GENIE v3.0.6
  // double BINDING_ENERGY = 0.0295; // 40Ar, GeV
  // This value is the shell-occupancy-weighted mean of the $E_{\alpha}$ values
  // listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an identical
  // procedure for 12C to obtain the binding energy value of 27.13 MeV, which is
  // adopted in their STV analysis described in arXiv:1910.08658.
  double BINDING_ENERGY = 0.02478; // 40Ar, GeV

  //useful masses
  double MASS_PROTON = 0.93827208; //GeV
  double MASS_MUON = 0.10565837; //GeV
  double MASS_PION0 = 0.13497666; //GeV
  double MASS_PIONPM =0.13957000; //GeV
  double MASS_NEUTRON = 0.93956541; // GeV

  //Counters                                                                                               
  ////////////////////////////////////////
  static const int num_cuts = 5;
  int n_events[num_cuts] = {0};
  int n_pmissT = 0;
  int bad_final = 0;


}
#endif
