///////////////////////
//7/19/2022: plotting.C
//Author: Samantha Sword-Fehlberg
//This makes plots that compare the selected events from 3 (4) GENIE MC models as a function of a number of variables.
//We previousl were making area normalized to compare the shapes, but this updated code now makes cross section plots as well.
//
//[HOW TO RUN THE CODE]
//root -b plotting.C
//plotting s
//s.main()
//
//PNG images of the plots are saved to specified pathname
///////////////////////
#include "plotting.h"
#include "tools/constants.h"
using namespace Constants;

void plotting::main(){

  //Pull the constants from constants. Determines which file we open
  const char* FSI;
  const char* XSEC_BINNING;

  if(fsi == true){FSI = "fsi";} else { FSI = "no_fsi";};
  if(use_xsec_binning == true){ XSEC_BINNING = "xsec";} else { XSEC_BINNING = "no_xsec";};

  //Load in the histogram file
  std::vector<TFile*> files; //vector of all the files
  TFile* f1 = new TFile(Form("/uboone/data/users/sfehlber/CC2p/MC_Studies/GENIE/histograms/%s/hists_empirical_lwellyn_%s_%s.root",FSI,FSI,XSEC_BINNING));
  TFile* f2 = new TFile(Form("/uboone/data/users/sfehlber/CC2p/MC_Studies/GENIE/histograms/%s/hists_nieves_%s_%s.root",FSI,FSI,XSEC_BINNING));
  TFile* f3 = new TFile(Form("/uboone/data/users/sfehlber/CC2p/MC_Studies/GENIE/histograms/%s/hists_susav2_%s_%s.root",FSI,FSI,XSEC_BINNING));
  //TFile* f4=new TFile("/uboone/data/users/sfehlber/CC2p/MC_Studies/GENIE/histograms/%s/hists_GCF_CCQE_%s%s.root",FSI,FSI,XSEC_BINNING)); //don't know how to compare to this yet. 
  files.push_back(f1);
  files.push_back(f2);
  files.push_back(f3);
  //files.push_back(f4);
    
  const char* pathname = Form("/uboone/data/users/sfehlber/CC2p/MC_Studies/images/%s/%s",FSI,XSEC_BINNING);
  string path(pathname);
   
  //Grab the histograms
  //////////////////////
  for(int i = 0; i < num_cuts; i++){
    for(int j = 0; j < num_mc; j++){    
      for(int k = 0; k < num_variables; k++){
	h_physics[j][i][k] = (TH1D*)files[j]->Get(Form("h%s%s",physics_variables[k],cuts[i]));//mcc8

      }
      for(int k = 0; k < num_var; k++){
	h_muon[j][i][k] = (TH1D*)files[j]->Get(Form("h_muon%s%s",var[k],cuts[i])); //mcc8 
	h_recoil[j][i][k] = (TH1D*)files[j]->Get(Form("h_recoil%s%s",var[k],cuts[i])); //mcc8 
	h_leading[j][i][k] = (TH1D*)files[j]->Get(Form("h_leading%s%s",var[k],cuts[i])); //mcc8
      }
    }
  }

  //These were things from Raquel's original analysis.
  //Don't use these now: 7/25/2022
  /*TH1D* h_Enu = (TH1D*)f1->Get("h_Enu");
  TH1D* h_pmissT = (TH1D*)f1->Get("h_pmissT");
  TH1D* h_pMiss = (TH1D*)f1->Get("h_pMiss");
  TH1D* h_cos_gamma_lab = (TH1D*)f1->Get("h_cos_gamma_lab");
  TH1D* h_cos_gamma_cm_with_cut = (TH1D*)f1->Get("h_cos_gamma_cm_with_cut");
  TH2D* h_2D_p1vp2 = (TH2D*)f1->Get("h_2D_p1vp2");
  TH2D* h_2D_cos_gamma_cm_vs_p_missT = (TH2D*)f1->Get("h_2D_cos_gamma_cm_v_p_missT");*/

  std::cout<<"Finished Grabbing Histograms"<<std::endl;
  
  //Plotting Time!
  ////////////////
  for(int i = 0; i < num_cuts; i++){

    //Physics variables
    for(int j =0; j < num_variables; j++){
      cross_section_plot(h_physics[0][i][j],h_physics[1][i][j],h_physics[2][i][j],y_lim_physics_xsec[i][j],titles_physics[j],path.c_str(),Form("%s%s",cuts[i],physics_variables[j]));
      area_normalized_plot(h_physics[0][i][j],h_physics[1][i][j],h_physics[2][i][j],y_lim_physics[i][j],titles_physics[j],path.c_str(),Form("%s%s",cuts[i],physics_variables[j]));
    } //end loop over num_variables

    //particles
    for(int j = 0; j < num_var; j++){

      //Muon
      cross_section_plot(h_muon[0][i][j],h_muon[1][i][j],h_muon[2][i][j],muon_ylim_xsec[j],Form("Muon %s",titles_var[j]),path.c_str(),Form("%s_muon%s",cuts[i],var[j]));
      area_normalized_plot(h_muon[0][i][j],h_muon[1][i][j],h_muon[2][i][j],muon_ylim[j],Form("Muon %s",titles_var[j]),path.c_str(),Form("%s_muon%s",cuts[i],var[j])); //muon
      
      //Leading Proton
      cross_section_plot(h_leading[0][i][j],h_leading[1][i][j],h_leading[2][i][j],leading_ylim_xsec[j],Form("Leading %s",titles_var[j]),path.c_str(),Form("%s_leading%s",cuts[i],var[j]));
      area_normalized_plot(h_leading[0][i][j],h_leading[1][i][j],h_leading[2][i][j],leading_ylim[j],Form("Leading Proton %s",titles_var[j]),path.c_str(),Form("%s_leading%s",cuts[i],var[j])); //leadin

      //Recoil Proton
      cross_section_plot(h_recoil[0][i][j],h_recoil[1][i][j],h_recoil[2][i][j],recoil_ylim_xsec[j],Form("Recoil %s",titles_var[j]),path.c_str(),Form("%s_recoil%s",cuts[i],var[j]));
      area_normalized_plot(h_recoil[0][i][j],h_recoil[1][i][j],h_recoil[2][i][j],recoil_ylim[j],Form("Recoil Proton %s",titles_var[j]),path.c_str(),Form("%s_recoil%s",cuts[i],var[j])); //recoil   

    } //end loop over number of var
  } //end loop over number of cuts

  //Now to make the Truth versus Pred. Plot for the Neutrino Energy and Pn,proxy
  /////////////////////////////////////////////////////////////////////////////

  for(int j = 0; j < num_mc; j++){ 
    h_true_nu_e[j] = (TH1D*)files[j]->Get(Form("h_nu_E_true%s",cuts[num_cuts-1]));
    h_true_pn[j] = (TH1D*)files[j]->Get(Form("h_pn_true%s",cuts[num_cuts-1])); 
  }

  truth_vs_pred_plot(h_physics[0][num_cuts-1][11],h_true_nu_e[0],
		     h_physics[1][num_cuts-1][11],h_true_nu_e[1],
		     h_physics[2][num_cuts-1][11],h_true_nu_e[2],
		     "#nu_{E} = E_{#mu} + T_{Lead} + T_{Recoil} + T_{A-2} + (P^{T}_{miss})^{2}/(2*M_{A-2})",
		     y_lim_physics[num_cuts-1][11],titles_physics[11],
		     path.c_str(),Form("%s%s",cuts[num_cuts-1],physics_variables[11])); //energy of the neutrino

  truth_vs_pred_plot(h_physics[0][num_cuts-1][7],h_true_pn[0],
                     h_physics[1][num_cuts-1][7],h_true_pn[1],
                     h_physics[2][num_cuts-1][7],h_true_pn[2],
                     "#p_{n,proxy} = #sqrt{p_{T}^{2} + p_{L}^{2}}",y_lim_physics[num_cuts-1][7],titles_physics[7],
                     path.c_str(),Form("%s%s",cuts[num_cuts-1],physics_variables[7])); //energy of the neutrino  

  //truth_vs_pred_plot(); //pn,proxy


  std::cout<<"----PROGRAM HAS FINISHED-----"<<std::endl;
}

