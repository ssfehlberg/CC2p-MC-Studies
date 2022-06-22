#include "helper_funcs.h"
using namespace Constants;
class variables{

 public:
  virtual vector<double> Calculate_Momenta(TVector3 vmuon, TVector3 vlead, TVector3 vrec);
  virtual vector<double> Calculate_Energy(TVector3 vmuon, TVector3 vlead, TVector3 vrec);
  virtual vector<double> Calculate_Detector_Angles(TVector3 vmuon, TVector3 vlead, TVector3 vrec);
  virtual vector<double> Calculate_Opening_Angles(TVector3 vmuon, TVector3 vlead, TVector3 vrec);
  virtual double Boost(TVector3 vLead, double ELead, TVector3 vRec, double ERec, TVector3 vMuon,double EMuon); //energies are total
  virtual vector<double> Calculate_STVS(bool add_protons, TVector3 vMuon, TVector3 vLead, TVector3 vRec);
  virtual double Calculate_Beam(TVector3 vMuon, TVector3 vLead, TVector3 vRec,double EMuon, double ELead, double ERec); //energies are kinetic 
  virtual void Calculate_Variables(TVector3 vmuon, TVector3 vlead, TVector3 vrec, bool add_protons);

  helper_funcs help; //class of helper functions

  vector<double> momenta;
  vector<double> Energies;
  vector<double> detector_angles;
  vector<double> opening_angles;
  double opening_angle_protons_COM;
  vector<double> stvs;
  double calculated_nu_E;

}; //End of class definition

//Calculates the Momenta of the three particles given the TVector3. Returns vector with magnitude of muon, lead, and reco momentum
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> variables::Calculate_Momenta(TVector3 vmuon, TVector3 vlead, TVector3 vrec){

  //muon
  double muon_mom;

  if((vmuon[0] == -9999. && vmuon[1] == -9999. && vmuon[2] == -9999.) || vmuon.Mag() > 9999.){
    muon_mom = -9999.0;
    std::cout<<"[VARIABLES] WE GOT A GARBAGE MOMENTUM FOR MUON!"<<std::endl;
  } else {
    muon_mom = double(vmuon.Mag());
  }
    
  //lead proton
  double lead_mom;
  if((vlead[0] == -9999. && vlead[1] == -9999. && vlead[2] == -9999.) || vlead.Mag() > 9999.){
    std::cout<<"[VARIABLES] WE GOT A GARBAGE MOMENTUM FOR LEAD!"<<std::endl;
    lead_mom = -9999.0;
  } else {
    lead_mom = double(vlead.Mag());
  }

  //recoil proton
  double rec_mom;
  if((vrec[0] == -9999. && vrec[1] == -9999. && vrec[2] == -9999.) || vrec.Mag() > 9999.){
    std::cout<<"[VARIABLES] WE GOT A GARBAGE MOMENTUM FOR RECOIL!"<<std::endl;
    rec_mom = -9999.0;
  } else {
    rec_mom = double(vrec.Mag());
  }

  //vector<double> momenta{muon_mom,lead_mom,rec_mom};
  momenta.push_back(muon_mom);
  momenta.push_back(lead_mom);
  momenta.push_back(rec_mom);
  return momenta;

}

//Calculates the Energies of the three particles given the TVector3. Returns vector with KE and Tot_E of muon, lead, and reco momentum                                                                                                              
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
vector<double> variables::Calculate_Energy(TVector3 vmuon, TVector3 vlead, TVector3 vrec){

  //muon
  double KE_Muon;
  double TotE_Muon;
  if(vmuon[0] == -9999. && vmuon[1] == -9999. && vmuon[2] == -9999.){
    KE_Muon = -9999.0;
    TotE_Muon = -9999.0;
  } else {
    KE_Muon = help.Which_Energy(false,vmuon,MASS_MUON);
    TotE_Muon = help.Which_Energy(true,vmuon,MASS_MUON);
  }

  //lead proton    
  double KE_Lead;
  double TotE_Lead;
  if(vlead[0] == -9999. && vlead[1] == -9999. && vlead[2] == -9999.){
    KE_Lead = -9999.0;
    TotE_Lead = -9999.0;
  } else {
    KE_Lead = help.Which_Energy(false,vlead,MASS_PROTON);
    TotE_Lead = help.Which_Energy(true,vlead,MASS_PROTON);
  }

  //recoil proton  
  double KE_Rec;
  double TotE_Rec;
  if(vrec[0] == -9999. && vrec[1] == -9999. && vrec[2] == -9999.){
    KE_Rec = -9999.0;
    TotE_Rec = -9999.0;
  } else {
    KE_Rec = help.Which_Energy(false,vrec,MASS_PROTON);
    TotE_Rec = help.Which_Energy(true,vrec,MASS_PROTON);
  }

  //vector<double> Energies{KE_Muon, TotE_Muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec};
  Energies.push_back(KE_Muon);
  Energies.push_back(TotE_Muon);
  Energies.push_back(KE_Lead);
  Energies.push_back(TotE_Lead);
  Energies.push_back(KE_Rec);
  Energies.push_back(TotE_Rec);
  return Energies;

} 

//Calculates the detector angles of the three particles given the TVector3. Returns vector with theta & phi of muon, lead, and reco            
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
vector<double> variables::Calculate_Detector_Angles(TVector3 vmuon, TVector3 vlead, TVector3 vrec){

  //muon
  double muon_phi;
  double muon_theta;
  if(vmuon[0] == -9999. && vmuon[1] == -9999. && vmuon[2] == -9999.){
    muon_phi = -9999.0;
    muon_theta = -9999.0;
  } else {
    muon_phi = double(vmuon.Phi());
    muon_theta = double(std::cos(vmuon.Theta()));
  }

  //lead proton                                                                                                                                                                                                                                                                    
  double lead_phi;
  double lead_theta;
  if(vlead[0] == -9999. && vlead[1] == -9999. && vlead[2] == -9999.){
    lead_phi = -9999.0;
    lead_theta = -9999.0;
  } else {
    lead_phi = double(vlead.Phi());
    lead_theta = double(std::cos(vlead.Theta()));
  }

  //recoil proton   
  double recoil_phi;
  double recoil_theta;                                                                                                                                                                                                          
  if(vrec[0] == -9999. && vrec[1] == -9999. && vrec[2] == -9999.){
    recoil_phi = -9999.0;
    recoil_theta = -9999.0;
  } else {
    recoil_phi = double(vrec.Phi());
    recoil_theta= double(std::cos(vrec.Theta()));
  }

  //vector<double> detector_angles{muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi};
  detector_angles.push_back(muon_theta);
  detector_angles.push_back(muon_phi);
  detector_angles.push_back(lead_theta);
  detector_angles.push_back(lead_phi);
  detector_angles.push_back(recoil_theta);
  detector_angles.push_back(recoil_phi);
  return detector_angles;

}

//Calculates the opening angles between protons, muon and lead proton, and muon and both protons all in lab frame given input TVector3. Returns opening angle protons, opening angle muonb+lead, opening angle muon_both         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
vector<double> variables::Calculate_Opening_Angles(TVector3 vmuon, TVector3 vlead, TVector3 vrec){

  //opening angle between the protons in the lab frame
  double opening_angle_protons_lab;
  if((vlead[0] == -9999. && vlead[1] == -9999. && vlead[2] == -9999.) || (vrec[0] == -9999. && vrec[1] == -9999. && vrec[2] == -9999.)){
    opening_angle_protons_lab = -9999.0;
  } else{
    opening_angle_protons_lab = std::cos(vlead.Angle(vrec)); //((vlead[0]*vrec[0])+(vlead[1]*vrec[1])+(vlead[2]*vrec[2]))/(vlead.Mag()*vrec.Mag()); //vlead.Angle(vrec); 
  }

  //opening angle between the muon and leading proton
  double opening_angle_protons_mu_leading;
  if(vlead[0] == -9999. && vlead[1] == -9999. && vlead[2] == -9999.){
    opening_angle_protons_mu_leading = -9999.0;
  }else{
    opening_angle_protons_mu_leading = std::cos(vmuon.Angle(vlead)); //((vlead[0]*vMuon[0])+(vlead[1]*vMuon[1])+(vlead[2]*vMuon[2]))/(vlead.Mag()*vMuon.Mag());                                                                           
  }
  
  //opening angle between the muon and both protons
  double opening_angle_protons_mu_both;
  TVector3 vProton(vlead[0]+vrec[0],vlead[1]+vrec[1],vlead[2]+vrec[2]); 
  if((vlead[0] == -9999. && vlead[1] == -9999. && vlead[2] == -9999.) || (vrec[0] == -9999. && vrec[1] == -9999. && vrec[2] == -9999.)){
    opening_angle_protons_mu_both = -9999.0;
  } else{
    opening_angle_protons_mu_both = std::cos(vmuon.Angle(vProton)); //((vProton[0]*vMuon[0])+(vProton[1]*vMuon[1])+(vProton[2]*vMuon[2]))/(vProton.Mag()*vMuon.Mag());
  }

  //vector<double> opening_angles{opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both};
  opening_angles.push_back(opening_angle_protons_lab);
  opening_angles.push_back(opening_angle_protons_mu_leading);
  opening_angles.push_back(opening_angle_protons_mu_both);
  return opening_angles;

 }

//Calculates the opening angle between the protons in the COM frame. INPUT ENERGIES MUST BE TOTAL ENERGIES!
///////////////////////////////////////////////////////////////////////////////////////////////////////////
double variables::Boost(TVector3 vLead, double ELead, TVector3 vRec, double ERec, TVector3 vMuon, double EMuon){

  TLorentzVector LEAD(vLead[0],vLead[1],vLead[2],ELead);
  TLorentzVector REC(vRec[0],vRec[1],vRec[2],ERec);
  TLorentzVector betacm(vRec[0]+vLead[0]+vMuon[0],vRec[1]+vLead[1]+vMuon[1],vRec[2]+vLead[2]+vMuon[2],ERec+ELead+EMuon);
  TVector3 boost = betacm.BoostVector(); //the boost vector                                                                        
  LEAD.Boost(-boost); //boost leading proton                                                                                                                                                                                                 
  REC.Boost(-boost); //boost recoil proton
  opening_angle_protons_COM = cos(LEAD.Angle(REC.Vect()));

  //Sanity check that I am doing the boost right
  if(_debug){
  
    TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);
    TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec);
    TLorentzVector muon(vMuon[0], vMuon[1], vMuon[2],EMuon);
    
    std::cout<<"[VARIABLES::BOOST] Value of all vectors added x: "<<lead[0]+rec[0]+muon[0]<<std::endl;
    std::cout<<"[VARIABLES::BOOST] Value of all vectors added y: "<<lead[1]+rec[1]+muon[1]<<std::endl;
    std::cout<<"[VARIABLES::BOOST] Value of all vectors added z: "<<lead[2]+rec[2]+muon[2]<<std::endl;
  
    TLorentzVector betacm_check(vRec[0]+vLead[0],vRec[1]+vLead[1],vRec[2]+vLead[2],ERec+ELead);
    TVector3 boost_check = betacm_check.BoostVector(); //the boost vector                                                 
    lead.Boost(-boost_check); //boost leading proton                                                                      
    rec.Boost(-boost_check); //boost recoil proton                                                                               
    std::cout<<"[VARIABLES::BOOST] Value of lead+rec boosted vectors x: "<<lead[0]+rec[0]<<std::endl;
    std::cout<<"[VARIABLES::BOOST] Value of lead+rec boosted vectors y: "<<lead[1]+rec[1]<<std::endl;
    std::cout<<"[VARIABLES::BOOST] Value of lead+rec boosted vectors z: "<<lead[2]+rec[2]<<std::endl;

  }

  return opening_angle_protons_COM;

}

//Calculates the STVs with given TVector3. add_protons indicates if we are adding the lead and rec proton momenta together.  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
vector<double> variables::Calculate_STVS(bool add_protons, TVector3 vMuon, TVector3 vLead, TVector3 vRec){

  //add the protons together (if you want of courses)
  TVector3 vProton;
  if(add_protons){
    vProton.SetXYZ(vLead[0]+vRec[0],vLead[1]+vRec[1],vLead[2]+vRec[2]);
  }else{
    vProton.SetXYZ(vLead[0],vLead[1],vLead[2]);
  }

  //delta PT and delta phiT
  double delta_pT;
  double delta_phiT;
  if((vLead[0] == -9999. && vLead[1] == -9999. && vLead[2] == -9999.) || (vRec[0] == -9999. && vRec[1] == -9999. && vRec[2] == -9999.)){
    delta_pT = -9999.0;
    delta_phiT = -9999.0;
  }else {
    delta_pT = double((vMuon + vProton).Perp()); //perp takes the magnitude as well. stv delta pt                                                                                                                              
    delta_phiT = double(180.0/3.14)*double(std::acos( (-vMuon.X()*vProton.X() - vMuon.Y()*vProton.Y()) / (vMuon.XYvector().Mod() * vProton.XYvector().Mod()))); //stv delta phi t  
  }

  //delta alpha T
  double delta_alphaT;
  TVector3 delta_pT_vec;
  delta_pT_vec.SetXYZ((vMuon[0]+vProton[0]),(vMuon[1]+vProton[1]),0.0);

  if((vLead[0] == -9999. && vLead[1] == -9999. && vLead[2] == -9999.) || (vRec[0] == -9999. && vRec[1] == -9999. && vRec[2] == -9999.)){
    delta_alphaT = -9999.0;
  } else {
    delta_alphaT = double(180.0/3.14)*double(std::acos( (-vMuon.X()*delta_pT_vec.X()- vMuon.Y()*delta_pT_vec.Y()) / (vMuon.XYvector().Mod() * delta_pT_vec.XYvector().Mod()) )); //stv delta alpha T    
  }

  //pn
  double delta_pL;
  double pn;
  float Emu = std::sqrt(std::pow(MASS_MUON, 2) + vMuon.Mag2());
  float Ep = std::sqrt(std::pow(MASS_PROTON, 2) + vProton.Mag2());
  float R = MASS_TARGET + vMuon.Z() + vProton.Z() - Emu - Ep;

  // Estimated mass of the final remnant nucleus (CCQE assumption)
  float mf = MASS_TARGET - (MASS_NEUTRON + MASS_PROTON) + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);
  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );

  //vector<double> stvs{delta_pT,delta_alphaT,delta_phiT};
  stvs.push_back(delta_pT);
  stvs.push_back(delta_alphaT);
  stvs.push_back(delta_phiT);
  stvs.push_back(pn);
  return stvs;

}

//Calculates the Neutrino Energy using Raquel's Formula. INPUT ENERGIES MUST BE KE!
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
double variables::Calculate_Beam(TVector3 vMuon, TVector3 vLead, TVector3 vRec,double EMuon, double ELead, double ERec){

  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  if((vLead[0] == -9999. && vLead[1] == -9999. && vLead[2] == -9999.) || (vRec[0] == -9999. && vRec[1] == -9999. && vRec[2] == -9999.)){
    calculated_nu_E = -9999.0;
  } else {
    //based on formula in Raquel's technote (Eq. 1) and the argoneut paper  
    // E_nu = Total Energy Muon + kinetic Energy Proton + Kinetic Energy Proton + (P^{T}_{Miss})/2*M_{A-2} + E_{Miss}
    //P^{T}_{Miss}  = p_mu^T + p_{Lead}^T + p_{recoil}^T i.e. transverse components of all the vectors
    // E_{miss} is estimated to be 30.4 MeV. M_{A-2} is taken to be 35.37 GeV. 
    calculated_nu_E = (EMuon+MASS_MUON) + ELead + ERec +((PT_miss.Mag2())/(2.0*35.37)) + 0.0304;;
  }
  if(_debug) std::cout<<"[VARIABLES] Value of Calculated_Nu_E: "<<calculated_nu_E<<std::endl;
  return calculated_nu_E;

} 

//Calculates all the variables of interest using all the above functions cause I couldn't think of a better method.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void variables::Calculate_Variables(TVector3 vmuon, TVector3 vlead, TVector3 vrec, bool add_protons){

  Calculate_Momenta(vmuon,vlead,vrec); //returns momenta(muon_mom,lead_mom,rec_mom);
  std::vector<double> Energies = Calculate_Energy(vmuon,vlead,vrec); ///returns Energies(KE_muon, TotE_muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec);
  Calculate_Detector_Angles(vmuon,vlead,vrec); //returns detector_angles(muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi);
  Calculate_Opening_Angles(vmuon,vlead,vrec); //returns opening_angles(opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both);
  Boost(vlead,Energies[3],vrec,Energies[5],vmuon,Energies[1]); // returns opening_angle_protons_COM. input energies are total;
  Calculate_STVS(add_protons, vmuon, vlead, vrec); //returns STVS(delta_pT,delta_alphaT,delta_phiT);
  Calculate_Beam(vmuon,vlead,vrec,Energies[0],Energies[2],Energies[4]);//returns calculated_nu_E. input energies are KE;

} 
