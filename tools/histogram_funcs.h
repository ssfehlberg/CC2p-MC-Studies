#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <fstream> 
#include "vector"
#include "variables.h"
#include "constants.h"
using namespace Constants;

class histogram_funcs
{

 public:
  virtual void Define_Histograms();
  //virtual void Define_Parameters(TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vLead1,TVector3 vRec, double wght);
  //virtual void compute_raquel(TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, double wght);
  virtual void Fill_Histograms(int i,TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec,  double nu_true,TVector3 pn_vec,double wght);
  virtual void Scale_Histograms();
  virtual void Write_Histograms(); //writes histograms. works for all samples


 private:
  
  //Histograms 
  /////////////////
  static const int num_cuts = 5;//7
  const char* cuts[num_cuts] = {"_b4_cuts","_pmiss_cut","_muon_cut","_lead_cut","_rec_cut"};//,"_pionpm_cut","_pion0_cut"};
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  int num_bins[num_var] = {40,50,30,10}; //50 first bin, 10 last bin                                                                                                                                                                                         
  double xlim_low_muon[num_var] = {0.1,0,-1.5,-3.15}; //0.2 normally first -1.5                                                                                                                                                                              
  double xlim_low_proton[num_var] = {0.3,0,-1.5,-3.15};
  double xlim_high_recoil[num_var] = {1.0,0.35,1.5,3.15};
  double xlim_high_leading[num_var] = {1.0,0.6,1.5,3.15}; //1.5 normally in first, 1.2                                                                                                                                                                       
  double xlim_high_muon[num_var]={1.2,1,1.5,3.15};
  const char* xlabel[num_var] ={"P [GeV/c]","E [GeV]","cos(#theta)","#phi [Rad]"};

  TH1D* h_muon[num_cuts][num_var];
  TH1D* h_recoil[num_cuts][num_var];
  TH1D* h_leading[num_cuts][num_var];
  TH1D* h_opening_angle_protons_lab[num_cuts]; //opening angle between the protons in lab frame
  TH1D* h_opening_angle_protons_com[num_cuts]; //opening angle between the protons in com frame
  TH1D* h_opening_angle_mu_leading[num_cuts]; //opening angle between the muon and lead proton
  TH1D* h_opening_angle_mu_both[num_cuts]; //opening angle between both protons and the muon
  TH1D* h_delta_PT[num_cuts]; //stvs delta PT
  TH1D* h_delta_alphaT[num_cuts]; //stvs delta alphaT
  TH1D* h_delta_phiT[num_cuts]; //stvs delta phiT

  TH1D* h_pn[num_cuts]; //neutron momentum estimate
  TH1D* h_pn_true[num_cuts];

  TH1D* h_nu_E[num_cuts]; //neutrino energy estimate
  TH1D* h_nu_E_true[num_cuts];

  TH1D* h_mom_struck_nuc[num_cuts]; //estimate for neutron momentum: NOT GOOD
  TH1D* h_tot_pz[num_cuts]; //total pz of the system
  TH1D* h_tot_E[num_cuts]; //total energy of the system
  TH1D* h_tot_E_minus_beam[num_cuts]; //the total energy minus the beam
  TH1D* h_PT_squared[num_cuts]; //the PT squared quantitiy in Raquel's note

  //vertex location
  TH1D* h_vtxx;
  TH1D* h_vtxy;
  TH1D* h_vtxz;
  TH1D* h_vtxt;

  //stuff from Raquels technote
  TH1D* h_pmissT;
  TH1D* h_cos_gamma_lab;
  TH1D* h_cos_gamma_cm_with_cut;  
  TH2D* h_2D_cos_gamma_cm_v_p_missT;
  TH2D* h_2D_p1vp2;

  //List of histograms to make writing easier
  vector<TH1*> h_list;
  vector<TH2D*> h_list_2D;

  //Variables Class
  variables variables; //variables class

  /*//Parameters:
   /////////////////
   double EMuon; // Muon Energy
   double ELead; // Leading Energy
   double ERec; // Recoil Energy
   double omega; // Energy transfer
   double QSq; // Q-Squared
   double muon_theta; //note this is cos(theta)
   double lead_theta; //note this is cos(theta)
   double rec_theta; //note this is cos(theta)
   double muon_phi; //this is actual angle
   double lead_phi; //this is actual angle
   double rec_phi; //this is actual angle
   double open_angle; //note this is the cos(opening angle)
   double open_angle_mu; //note this is the cos(opening angle)  
   double open_angle_mu_proton;
   double delta_pT; //stv delta_pT
   double delta_alphaT; //stv delta_alphaT
   double delta_phiT; //stv delta_phiT
   double delta_phiT_sam; //to make the plot look better
   double delta_pL; //stv delta pL
   double pn; //reminant momentum ....?
   double p_missT; //transverse missing momentum
   double cos_gamma_lab; //cos(opening angle) in lab
   double cos_gamma_cm; //cos(opening angle) in cm
   double p_struck_nuc;
   double Enuc; //energy of struck nucleon
   double betacm; //boost parameter
   double pz_tot;
   double E_tot;
   double E_tot_minus_beam;
   double Eneutrino;
   TVector3 PT_miss;
   TVector3 vq;
   TVector3 vmiss;
   double Ev;
  */
  
}; //end of class


void histogram_funcs::Define_Histograms(){


  //vertex location
  h_vtxx = new TH1D("h_vtxx","h_vtxx; X Location of Vertex (cm); Counts",40,0,275);
  h_vtxy =new TH1D("h_vtxy","h_vtxy; Y Location of Vertex (cm); Counts",40,-125,125);
  h_vtxz = new TH1D("h_vtxz","h_vtxz; Z Location of Vertex (cm); Counts",50,0,1050);
  h_vtxt = new TH1D("h_vtxt","h_vtxt; Time of Interaction (#mus); Counts",50,2,10);
  h_list.push_back(h_vtxx);
  h_list.push_back(h_vtxy);
  h_list.push_back(h_vtxz);
  h_list.push_back(h_vtxt);
  
  //raquels stuff
  h_pmissT = new TH1D("h_pmissT","h_pmissT; P_{miss}^{T} (GeV/C)",50,0,2);
  h_cos_gamma_lab = new TH1D("h_cos_gamma_lab","h_cos_gamma_lab; cos(#gamma_{Lab}); Counts",30,-1.5,1.5);
  h_cos_gamma_cm_with_cut = new TH1D("h_cos_gamma_cm_with_cut","h_cos_gamma_cm_with_cut; cos(#gamma_{COM}); Counts",30,-1.5,1.5);
  h_2D_p1vp2 = new TH2D("h_2D_p1vp2","h_2D_p1vp2; Momentum of P_{Recoil} (Gev/c); Momentum of P_{Leading} (Gev/c)",50,0,3,50,0,3);
  h_2D_cos_gamma_cm_v_p_missT = new TH2D("h_2D_cos_gamma_cm_v_p_missT","h_2D_cos_gamma_cm_v_p_missT; P^{T}_{miss} (GeV/c); cos(#gamma_{COM})",20,0,1.6,20,-1.0,1.0);
  h_list.push_back(h_pmissT);
  h_list.push_back(h_cos_gamma_lab);
  h_list.push_back(h_cos_gamma_cm_with_cut);
  h_list.push_back(h_vtxx);
  h_list.push_back(h_vtxy);
  h_list.push_back(h_vtxz);
  h_list.push_back(h_vtxt);
  h_list_2D.push_back(h_2D_p1vp2);
  h_list_2D.push_back(h_2D_cos_gamma_cm_v_p_missT);

  for(int i = 0; i < num_cuts; i++){
    for(int j = 0; j < num_var; j++){
      if(use_xsec_binning == true){
	if(j == 0){ //momentum
	  //muon
	  const Int_t bins_muon = 5;
	  Double_t edges_muon[bins_muon+1] = {0.1,0.3,0.5,0.7,0.9,1.2};
	  h_muon[i][j] = new TH1D(Form("h_muon%s%s",var[j],cuts[i]),Form(" h_muon%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),bins_muon,edges_muon);
	  //recoil
	  const Int_t bins_recoil = 6;
	  Double_t edges_recoil[bins_recoil+1] = {0.3,0.4,0.5,0.6,0.7,0.8,1.0};
	  h_recoil[i][j] = new TH1D(Form("h_recoil%s%s",var[j],cuts[i]),Form("h_recoil%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),bins_recoil,edges_recoil);
	  //leading
	  const Int_t bins_leading = 6;
	  Double_t edges_leading[bins_leading+1] = {0.3,0.4,0.5,0.6,0.7,0.8,1.0}; 
	  h_leading[i][j] = new TH1D(Form("h_leading%s%s",var[j],cuts[i]),Form("h_leading%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),bins_leading,edges_leading);
	  
	}else if (j == 2){ //theta
	  const Int_t bins_theta = 8;
	  Double_t edges_theta[bins_theta+1] = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0};
	  h_muon[i][j] = new TH1D(Form("h_muon%s%s",var[j],cuts[i]),Form(" h_muon%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),bins_theta,edges_theta);
	  h_recoil[i][j] = new TH1D(Form("h_recoil%s%s",var[j],cuts[i]),Form("h_recoil%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),bins_theta,edges_theta);
	  h_leading[i][j] = new TH1D(Form("h_leading%s%s",var[j],cuts[i]),Form("h_leading%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),bins_theta,edges_theta);

	}else{ //phi and energy
	  h_muon[i][j] = new TH1D(Form("h_muon%s%s",var[j],cuts[i]),Form(" h_muon%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
	  h_recoil[i][j] = new TH1D(Form("h_recoil%s%s",var[j],cuts[i]),Form("h_recoil%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
	  h_leading[i][j] = new TH1D(Form("h_leading%s%s",var[j],cuts[i]),Form("h_leading%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
	} //end of phi
      
      }else if (use_xsec_binning == false){
	h_muon[i][j] = new TH1D(Form("h_muon%s%s",var[j],cuts[i]),Form(" h_muon%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),num_bins[j],xlim_low_muon[j],xlim_high_muon[j]);
	h_recoil[i][j] = new TH1D(Form("h_recoil%s%s",var[j],cuts[i]),Form("h_recoil%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),num_bins[j],xlim_low_proton[j],xlim_high_recoil[j]);
	h_leading[i][j] = new TH1D(Form("h_leading%s%s",var[j],cuts[i]),Form("h_leading%s%s ;%s; Counts",var[j],xlabel[j],cuts[i]),num_bins[j],xlim_low_proton[j],xlim_high_leading[j]);
      } //end of false statement
      
      h_list.push_back(h_muon[i][j]);
      h_list.push_back(h_recoil[i][j]);
      h_list.push_back(h_leading[i][j]);
    } //end of loop over variables for particles


    if(use_xsec_binning == true){
      const Int_t bins = 8;
      Double_t edges[bins+1] = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0};
      h_opening_angle_protons_lab[i] = new TH1D(Form("h_opening_angle_protons_lab%s",cuts[i]),Form("h_opening_angle_protons_lab%s; Opening Angle btwn Two Protons; Counts",cuts[i]),bins,edges);  
      h_opening_angle_protons_com[i]  = new TH1D(Form("h_opening_angle_protons_com%s",cuts[i]),Form("h_opening_angle_protons_com%s;cos(#gamma_{COM});Counts",cuts[i]),bins,edges);
      h_opening_angle_mu_leading[i]  = new TH1D(Form("h_opening_angle_mu_leading%s",cuts[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",cuts[i]),bins,edges);
      h_opening_angle_mu_both[i]  = new TH1D(Form("h_opening_angle_mu_both%s",cuts[i]),Form("h_opening_angle_mu_both%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",cuts[i]),bins,edges);

      const Int_t bins_stv_mom =7;
      Double_t edges_stv_mom[bins_stv_mom+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0};
      h_delta_PT[i]  = new TH1D(Form("h_delta_PT%s",cuts[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",cuts[i]),bins_stv_mom,edges_stv_mom);
      
      const Int_t bins_stv_angles = 6;
      Double_t edges_stv_angles[bins_stv_angles+1] = {0,30,60,90,120,150,180};
      h_delta_alphaT[i]  = new TH1D(Form("h_delta_alphaT%s",cuts[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",cuts[i]),bins_stv_angles,edges_stv_angles); 
      h_delta_phiT[i]  = new TH1D(Form("h_delta_phiT%s",cuts[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",cuts[i]),bins_stv_angles,edges_stv_angles);

      const Int_t bins_neutron = 8;
      Double_t edges_neutron[bins_neutron+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};
      h_pn[i]  = new TH1D(Form("h_pn%s",cuts[i]),Form("h_pn%s; p_{n} [GeV/c];Counts",cuts[i]),bins_neutron,edges_neutron);
      h_pn_true[i]  = new TH1D(Form("h_pn_true%s",cuts[i]),Form("h_pn_true%s; p_{n} [GeV/c];Counts",cuts[i]),bins_neutron,edges_neutron);

      const Int_t bins_nuE = 4;
      Double_t edges_nuE[bins_nuE+1] = {0.0, 0.4, 0.8, 1.2, 1.8};
      h_nu_E[i]  = new TH1D(Form("h_nu_E%s",cuts[i]),Form("h_nu_E%s; Total Energy; Counts;",cuts[i]),bins_nuE,edges_nuE);
      h_nu_E_true[i]  = new TH1D(Form("h_nu_E_true%s",cuts[i]),Form("h_nu_E_true%s; Total Energy; Counts;",cuts[i]),bins_nuE,edges_nuE);

      //basic bitches
      h_mom_struck_nuc[i]  = new TH1D(Form("h_mom_struck_nuc%s",cuts[i]),Form("h_mom_struck_nuc%s; P_{Init}; Counts", cuts[i]),30, 0, 1);
      h_tot_pz[i]  = new TH1D(Form("h_tot_pz%s",cuts[i]),Form("h_tot_pz%s; P_{Z}^{Total}; Counts",cuts[i]), 20, 0, 2);
      h_tot_E[i]  = new TH1D(Form("h_tot_E%s",cuts[i]),Form("h_tot_E%s; Total Energy; Counts;",cuts[i]),50,0,2.5);
      h_tot_E_minus_beam[i]  = new TH1D(Form("h_tot_E_minus_beam%s",cuts[i]),Form("h_tot_E_minus_beam%s; Total Energy Remaining (MeV/c); Counts;",cuts[i]),100,-100,0);
      h_PT_squared[i]  = new TH1D(Form("h_PT_squared%s",cuts[i]),Form("h_PT_squared%s; P_{T}^{2}; Counts", cuts[i]),50,0,5);
    
    } else if (use_xsec_binning == false){
      h_opening_angle_protons_lab[i]  = new TH1D(Form("h_opening_angle_protons_lab%s",cuts[i]),Form("h_opening_angle_protons_lab%s; Opening Angle btwn Two Protons; Counts",cuts[i]),30,-1.5,1.5); //50, 0, 1.5                       
      h_opening_angle_protons_com[i]  = new TH1D(Form("h_opening_angle_protons_com%s",cuts[i]),Form("h_opening_angle_protons_com%s;cos(#gamma_{COM});Counts",cuts[i]),30,-1.5,1.5);
      h_opening_angle_mu_leading[i]  = new TH1D(Form("h_opening_angle_mu_leading%s",cuts[i]),Form("h_opening_angle_mu_leading%s;Opening Angle btwn Muon and Leading Proton; Counts",cuts[i]),30,-1.5,1.5);
      h_opening_angle_mu_both[i]  = new TH1D(Form("h_opening_angle_mu_both%s",cuts[i]),Form("h_opening_angle_mu_both%s; Opening Angle btwn Muon and Total Proton Momentum; Counts",cuts[i]),30,-1.5,1.5);
      h_delta_PT[i]  = new TH1D(Form("h_delta_PT%s",cuts[i]),Form("h_deltaPT%s;#delta P_{T} [GeV/c];Counts",cuts[i]),15,0,1); //normally 10 bins                                                                                      
      h_delta_alphaT[i]  = new TH1D(Form("h_delta_alphaT%s",cuts[i]),Form("h_delta_alphaT%s; #delta #alpha_{T} [Deg.];Counts",cuts[i]),10,0,180); //0,180                                                                            
      h_delta_phiT[i]  = new TH1D(Form("h_delta_phiT%s",cuts[i]),Form("h_delta_phiT%s; #delta #phi_{T} [Deg.];Counts",cuts[i]),10,0,180); //0,180                                                                                     
      h_pn[i]  = new TH1D(Form("h_pn%s",cuts[i]),Form("h_pn%s; p_{n} [GeV/c];Counts",cuts[i]),25,0,0.5);
      h_pn_true[i]  = new TH1D(Form("h_pn_true%s",cuts[i]),Form("h_pn_true%s; p_{n} [GeV/c];Counts",cuts[i]),25,0,0.5);

      h_nu_E[i]  = new TH1D(Form("h_nu_E%s",cuts[i]),Form("h_nu_E%s; Total Energy; Counts;",cuts[i]),50,0,2.5);
      h_nu_E_true[i]  = new TH1D(Form("h_nu_E_true%s",cuts[i]),Form("h_nu_E_true%s; Total Energy; Counts;",cuts[i]),50,0,2.5);

      h_mom_struck_nuc[i]  = new TH1D(Form("h_mom_struck_nuc%s",cuts[i]),Form("h_mom_struck_nuc%s; P_{Init}; Counts", cuts[i]),30, 0, 1);
      h_tot_pz[i]  = new TH1D(Form("h_tot_pz%s",cuts[i]),Form("h_tot_pz%s; P_{Z}^{Total}; Counts",cuts[i]), 20, 0, 2);
      h_tot_E[i]  = new TH1D(Form("h_tot_E%s",cuts[i]),Form("h_tot_E%s; Total Energy; Counts;",cuts[i]),50,0,2.5);
      h_tot_E_minus_beam[i]  = new TH1D(Form("h_tot_E_minus_beam%s",cuts[i]),Form("h_tot_E_minus_beam%s; Total Energy Remaining (MeV/c); Counts;",cuts[i]),100,-100,0);
      h_PT_squared[i]  = new TH1D(Form("h_PT_squared%s",cuts[i]),Form("h_PT_squared%s; P_{T}^{2}; Counts", cuts[i]),50,0,5);
    } //end loop over false statemment

    h_list.push_back(h_PT_squared[i] );
    h_list.push_back(h_opening_angle_mu_both[i] );
    h_list.push_back(h_nu_E[i]);
    h_list.push_back(h_nu_E_true[i]);
    h_list.push_back(h_tot_E[i]);
    h_list.push_back(h_tot_E_minus_beam[i] );
    h_list.push_back(h_opening_angle_protons_lab[i] );
    h_list.push_back(h_opening_angle_protons_com[i] );
    h_list.push_back(h_opening_angle_mu_leading[i] );
    h_list.push_back(h_delta_PT[i] );
    h_list.push_back(h_delta_alphaT[i] );
    h_list.push_back(h_delta_phiT[i] );
    h_list.push_back(h_pn[i] );
    h_list.push_back(h_pn_true[i] );
    h_list.push_back(h_mom_struck_nuc[i] );
    h_list.push_back(h_tot_pz[i] );

  }//end loop over number of cuts

  // Set each histogram to "Sumw2" so that weights are handled correctly
  for (int i = 0; i < h_list.size(); i++){
    h_list[i]->Sumw2();
  }
  for(int i = 0; i < h_list_2D.size(); i++){
    h_list_2D[i]->Sumw2();
  }
} //end of define histograms

/*
void histogram_funcs::Define_Parameters(TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vLead1, TVector3 vRec, double wght){ 
  //General stuff
  EMuon = sqrt(MUON_MASS*MUON_MASS + vMuon.Mag2()) - MUON_MASS; // Muon Energy
  ELead = sqrt(PROTON_MASS*PROTON_MASS + vLead.Mag2()) - PROTON_MASS; // Leading Energy
  ERec = sqrt(PROTON_MASS*PROTON_MASS + vRec.Mag2()) - PROTON_MASS; // Recoil Energy
  E_tot = (EMuon + MUON_MASS) + ELead + ERec;

  //Beam Stuff
  PT_miss.SetXYZ(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  Eneutrino = (EMuon+MUON_MASS) + ELead + ERec +((PT_miss.Mag2())/(2.0*35.37)) + 0.0304;
  //TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                            
  vq = vBeam - vMuon; // Momentum transfer                                                                  
  vmiss = vLead - vq; // Missing momentum        
  E_tot_minus_beam = (E_tot - Eneutrino) * 1000;
  TVector3 vProton;
  if(add_protons){
    vProton.SetXYZ(vLead[0]+vRec[0],vLead[1]+vRec[1],vLead[2]+vRec[2]);
  }else{
    vProton.SetXYZ(vLead[0],vLead[1],vLead[2]);
  }

  Ev = vBeam[2];
  omega = Ev - EMuon; // Energy transfer
  QSq = vq.Mag2() - omega*omega; // Q-Squared
  muon_theta = vMuon[2]/vMuon.Mag(); //note this is cos(theta)
  lead_theta = vLead[2]/vLead.Mag(); //note this is cos(theta)
  rec_theta = vRec[2]/vRec.Mag(); //note this is cos(theta)
  muon_phi = atan2(vMuon[0],vMuon[1]); //this is actual angle
  lead_phi =atan2(vLead[0],vLead[1]); //this is actual angle
  rec_phi =atan2(vRec[0],vRec[1]); //this i actual angle
  open_angle = ((vLead[0]*vRec[0])+(vLead[1]*vRec[1])+(vLead[2]*vRec[2]))/(vLead.Mag()*vRec.Mag()); //note this is the cos(opening angle)   
  open_angle_mu = ((vLead[0]*vMuon[0])+(vLead[1]*vMuon[1])+(vLead[2]*vMuon[2]))/(vLead.Mag()*vMuon.Mag()); //note this is the cos(opening angle)
   open_angle_mu_proton = ((vProton[0]*vMuon[0])+(vProton[1]*vMuon[1])+(vProton[2]*vMuon[2]))/(vProton.Mag()*vMuon.Mag()); //cos(opening angle) between the total proton momentum vector and the muon


  //Stuff for STVs
  delta_pT = (vMuon + vLead1).Perp();
  delta_phiT = std::acos( (-vMuon.X()*vLead1.X() - vMuon.Y()*vLead1.Y()) / (vMuon.XYvector().Mod() * vLead1.XYvector().Mod()));
  TVector2 delta_pT_vec = (vMuon + vLead1).XYvector();
  delta_alphaT = std::acos( (-vMuon.X()*delta_pT_vec.X()- vMuon.Y()*delta_pT_vec.Y()) / (vMuon.XYvector().Mod() * delta_pT_vec.Mod()) );
  double Emu = std::sqrt(std::pow(MUON_MASS, 2) + vMuon.Mag2());
  double Ep = std::sqrt(std::pow(PROTON_MASS, 2) + vLead.Mag2());
  double R = TARGET_MASS + vMuon.Z() + vLead.Z() - Emu - Ep;
  double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY; //Estimated mass of the final remnant nucleus (CCQE assumption))
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);
  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );
  
  //Stuff for Raquel
  p_missT = std::sqrt(std::pow((vMuon[0]+vLead[0]+vRec[0]),2)+std::pow(((vMuon[1]+vLead[1]+vRec[1])),2)); //transverse missing momentum 
  cos_gamma_lab = cos(vLead.Angle(vRec)); //cos(opening angle in lab frame) but fancy!
  Enuc = std::sqrt(std::pow(NEUTRON_MASS,2) + vmiss.Mag2()); //energy of struck nucleon
  TLorentzVector lead(vLead[0],vLead[1],vLead[2],ELead);//leading proton TLorentzVector
  TLorentzVector rec(vRec[0],vRec[1],vRec[2],ERec); //Recoil proton TLorentzVector
  TLorentzVector betacm(vmiss[0] + vRec[0] + vBeam[0],vmiss[1] + vRec[1] + vBeam[1], vmiss[2] + vRec[2]+ vBeam[2], Enuc + ERec + Ev); //beta for CM
  TVector3 boost = betacm.BoostVector(); //the boost vector
  lead.Boost(-boost); //boost leading proton
  rec.Boost(-boost); //boost recoil proton
  cos_gamma_cm = cos(lead.Angle(rec.Vect())); //uses Lorentz Vectors

  //Struck nucleon Momentum:
  TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
  TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
  p_struck_nuc = p_struck_nuc_vector.Mag();
  pz_tot = vLead[2] + vRec[2];

  //Fill da random histograms
  h_pMiss->Fill(vmiss.Mag(), wght); // Filling the histogram of missing mass magnitude
  h_QSq->Fill(QSq, wght); //the q2
  h_Enu->Fill(Ev,wght);

}

void histogram_funcs::compute_raquel(TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, double wght){
  h_2D_p1vp2->Fill(vRec.Mag(),vLead.Mag(),wght); //figure 1 in Raquels note
  h_pmissT->Fill(p_missT,wght); //plot of the transverse missing momentum for my sake
  h_cos_gamma_lab->Fill(cos_gamma_lab,wght); //cos(gamma of lab frame) for sanity check
  h_2D_cos_gamma_cm_v_p_missT->Fill(p_missT,cos_gamma_cm,wght); //Figure 3 in Raquel's Note  
  
  //Remove events with P_missT greater than 300 MeV
  //if(p_missT >= 0.3) continue; //this is broken!!!!!!!
  //n_pmissT++;
  h_cos_gamma_cm_with_cut->Fill(cos_gamma_cm,wght); //Figure 4 in Raquel's Note

  //Figure 5 in Raquel's Note: Dont' know how to get P perp pn-p2
  //Figure 7 in Raquel's Note //need cut on previous plot
}
*/
void histogram_funcs::Fill_Histograms(int i,TVector3 vBeam, TVector3 vMuon, TVector3 vLead, TVector3 vRec, double nu_true, TVector3 pn_vec, double wght){

    // Run the Calculate Variables function inside of variables.h. Returns following:
  // 1) vector: momenta(muon_mom,lead_mom,rec_mom);
  // 2) vector: Energies(KE_muon, TotE_muon, KE_Lead, TotE_Lead, KE_Rec, TotE_Rec);    
  // 3) vector: detector_angles(muon_theta,muon_phi,lead_theta,lead_phi,recoil_theta,recoil_phi);
  // 4) vector: opening_angles(opening_angle_protons_lab,opening_angle_protons_mu_leading,opening_angle_protons_mu_both);  // 5) double: opening_angle_protons_COM 
  // 6) vector: STVS(delta_pT,delta_alphaT,delta_phiT); 
  // 7) double: calculated_nu_E

  variables.Calculate_Variables(vMuon,vLead,vRec,add_protons);

  //muon
  double muon_mom = variables.momenta[0];
  double muon_theta = variables.detector_angles[0]; //cosine applied
  double muon_phi = variables.detector_angles[1];
  double EMuon = variables.Energies[0]; //KE
  
  //lead proton
  double lead_mom = variables.momenta[1];
  double lead_theta = variables.detector_angles[2]; //cosine applied
  double lead_phi = variables.detector_angles[3];
  double ELead = variables.Energies[2]; //KE

  //recoil proton
  double recoil_mom = variables.momenta[2];
  double recoil_theta = variables.detector_angles[4]; //cosine applied
  double recoil_phi =  variables.detector_angles[5];
  double ERec = variables.Energies[4]; //KE

  //opening angles
  double opening_angle_protons_lab = variables.opening_angles[0]; //cosine applied
  double opening_angle_protons_COM = variables.opening_angle_protons_COM; //cosine applied
  double opening_angle_protons_mu_leading = variables.opening_angles[1]; //cosine applied
  double opening_angle_protons_mu_both = variables.opening_angles[2]; //cosine applied

  //stvs
  double delta_PT = variables.stvs[0];
  double delta_alphaT = variables.stvs[1]; //degrees
  double delta_phiT = variables.stvs[2]; //degrees
  double pn = variables.stvs[3]; //GeV/c

  //neutrino energy & PT Miss
  TVector3 PT_miss(vMuon[0]+vLead[0]+vRec[0],vMuon[1]+vRec[1]+vLead[1],0);
  double Eneutrino = variables.calculated_nu_E;
  double Eneutrino_true = nu_true;

  //Beam Stuff
  double E_tot = (EMuon + MASS_MUON) + ELead + ERec;
  //TVector3 vBeam(0.,0.,Eneutrino); // z-direction is defined along the neutrino direction                            
  TVector3 vq = vBeam - vMuon; // Momentum transfer                                                                  
  TVector3 vmiss = vLead - vq; // Missing momentum        
  double E_tot_minus_beam = (E_tot - Eneutrino) * 1000;
  
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of PT_miss magnitude: "<<PT_miss.Mag()<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of PT_miss magnitude2: "<<PT_miss.Mag2()<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of PT_miss magnitude2 divided by 2*35.37: "<<(PT_miss.Mag2())/(2.0*35.37)<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of Eneutrino: "<<Eneutrino<<std::endl;
  if(_debug) std::cout<<"[HISTOGRAM_FUNCS] Value of E_tot_minus_beam: "<<E_tot_minus_beam<<std::endl;

  //Struck nucleon Momentum:                                                                                            
  TVector3 vector_sum(vMuon[0] + vLead[0] + vRec[0], vMuon[1] + vLead[1] + vRec[1], vMuon[2] + vLead[2] + vRec[2]);
  TVector3 p_struck_nuc_vector(vector_sum[0], vector_sum[1], 0);
  double p_struck_nuc = p_struck_nuc_vector.Mag();
  double pz_tot = vLead[2] + vRec[2];

  //Energy of Struck Nucleon:
  double En = std::sqrt(std::pow(MASS_NEUTRON,2) + vmiss.Mag2()); //energy of struck nucleon 

  h_muon[i][0]->Fill(muon_mom,wght);
  h_muon[i][1]->Fill(EMuon,wght);
  h_muon[i][2]->Fill(muon_theta,wght);
  h_muon[i][3]->Fill(muon_phi,wght);

  h_leading[i][0]->Fill(vLead.Mag(),wght);
  h_leading[i][1]->Fill(ELead,wght);
  h_leading[i][2]->Fill(lead_theta,wght);
  h_leading[i][3]->Fill(lead_phi,wght);

  h_recoil[i][0]->Fill(vRec.Mag(),wght);
  h_recoil[i][1]->Fill(ERec,wght);
  h_recoil[i][2]->Fill(recoil_theta,wght);
  h_recoil[i][3]->Fill(recoil_phi,wght);

  h_opening_angle_protons_lab[i]->Fill(opening_angle_protons_lab,wght);
  h_opening_angle_protons_com[i]->Fill(opening_angle_protons_COM,wght);
  h_opening_angle_mu_leading[i]->Fill(opening_angle_protons_mu_leading,wght);
  h_opening_angle_mu_both[i]->Fill(opening_angle_protons_mu_both,wght);
  h_delta_PT[i]->Fill(delta_PT,wght);
  h_delta_alphaT[i]->Fill(delta_alphaT,wght);
  h_delta_phiT[i]->Fill(delta_phiT,wght);

  //TVector3 pn_vec(pnx,pny,pny);
  h_pn[i]->Fill(pn,wght);
  h_pn_true[i]->Fill(pn_vec.Mag(),wght);

  h_nu_E[i]->Fill(Eneutrino,wght);
  h_nu_E_true[i]->Fill(Eneutrino_true,wght);

  h_mom_struck_nuc[i]->Fill(p_struck_nuc,wght);
  h_tot_pz[i]->Fill(pz_tot,wght);
  h_tot_E[i]->Fill(E_tot,wght);
  h_tot_E_minus_beam[i]->Fill(E_tot_minus_beam,wght);
  h_PT_squared[i]->Fill(PT_miss.Mag2(),wght);

  //have to make sure to clear the variables before you leave
  //////////////////////////////////////////////////////////
  variables.momenta.clear();
  variables.detector_angles.clear();
  variables.opening_angles.clear();
  variables.stvs.clear();
  variables.Energies.clear();

}

void histogram_funcs::Scale_Histograms(){
  int n_bins = 0;
  double n_events = 0;
  double bin_content = 0;
  double value = 0;

  for(int i=0 ; i < h_list.size(); i++){
    n_bins = h_list[i]->GetNbinsX();
    n_events = h_list[i]->Integral();
   
    for(int j=1; j < n_bins + 1; j++){
      bin_content = h_list[i]->GetBinContent(j);
      value = bin_content/n_events;
      h_list[i]->SetBinContent(j, value);
    }
  }

}

void histogram_funcs::Write_Histograms(){ 

  for(int i = 0; i < h_list.size(); i++){
    h_list[i]->Write();
  }
   for(int i = 0; i < h_list_2D.size(); i++){
     h_list_2D[i]->Write();
   }
  
}
