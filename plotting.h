/////////////////////
//7/19/2022: plotting.h
//Author: Samantha Sword-Fehlberg
//Header file for plotting.C.
//
////////////////////////////

//Helpful includes
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <string>

class plotting{

 public:
  virtual void main();
  virtual void area_normalized_plot(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, double y_lim,const char* title, const char* path, const char* variable);
  virtual void cross_section_plot(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF, double y_lim,const char* title, const char* path, const char* variable);
  virtual void truth_pred_plot();

 private:

  //Stuff for date and time
  //////////////////////////
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;

  //Generally stuff for all the plots
  //////////////////////////////////
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);
  const char* in_prog= "#scale[0.6]{MicroBooNE In-Progress}";
  const char* prelim= "#scale[0.6]{MicroBooNE Preliminary}";

  //Cut specifications
  ////////////////////////
  static const int num_mc = 3; //Number of MC samples. Normally 4, but don't use the GCF sample right now
  static const int num_cuts = 5;//Number of cuts applied. It is also possible to include cuts on the pion momentum (7), but I don't include them here (5). 
  const char* cuts[num_cuts] = {"_b4_cuts","_pmiss_cut","_muon_cut","_lead_cut","_rec_cut"};//names of the cuts in the root file
  //const char* cuts[num_cuts] = {"_b4_cuts","_pmiss_cut","_muon_cut","_lead_cut","_rec_cut","_pionpm_cut","_pion0_cut"};                                                                                                                                                                              
  const char* titles_cuts[num_cuts] = {"Before Cuts","After P_{Miss} Cut","After Muon Mom. Cut","After Leading Mom. Cut","After Recoil Mom. Cut"}; //Pretty version of the cut names
  
  //Individual particle plots
  ///////////////////////////
  static const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  const char* titles_var[num_var] = {"Momentum (GeV/c)","Energy (GeV)","cos(#theta)","#phi (Rad)"};
  TH1D* h_muon[num_mc][num_cuts][num_var];
  TH1D* h_leading[num_mc][num_cuts][num_var];
  TH1D* h_recoil[num_mc][num_cuts][num_var];
  double muon_ylim[num_var] = {0.09,0.12,0.6,0.05}; //y limit for area normalization plot       
  double muon_ylim_xsec[num_var] = {0.09,0.12,0.6,0.05}; //y limit for xsec plot     
  double leading_ylim[num_var] = {0.1,0.3,0.6,0.05}; //y limit for area normalization plot
  double leading_ylim_xsec[num_var] = {0.09,0.12,0.6,0.05}; //y limit for xsec plot
  double recoil_ylim[num_var] = {0.25,0.5,0.3,0.05}; //y limit for area normalization plot 
  double recoil_ylim_xsec[num_var] = {0.09,0.12,0.6,0.05}; //y limit for xsec plot

  //crap ton of physics plots
  /////////////////////////////
  static const int num_variables = 12;
  const char* physics_variables[num_variables] = {"_opening_angle_protons_lab","_opening_angle_protons_com",
						  "_opening_angle_mu_leading","_opening_angle_mu_both",
						  "_delta_PT","_delta_alphaT","_delta_phiT","_pn",
						  "_tot_E","_tot_E_minus_beam","_PT_squared","_nu_E"};//_true_minus_calculated_E_neutrino"};
  const char* titles_physics[num_variables] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})",
					       "cos(#gamma_{#mu, p_{L}})", "cos(#gamma_{#mu, p_{L}+P_{R}})",
					       "#deltaP_{T} (GeV/c)","#delta#alpha_{T} (Deg.)","#delta#phi_{T} (Deg.)","p_{n} (GeV/c)",
					       "Total Muon Energy + Total Kinetic Energy of Protons (GeV/c)","Total Energy Minus Beam Energy (MeV/c)","(P^{T}_{miss})^{2} (GeV^{2}/c^{2})","Energy of the Neutrino (GeV/c)"};//,"E^{True}_{#nu} - E^{Calculated}_{#nu} (GeV/c)"};
  TH1D* h_physics[num_mc][num_cuts][num_variables];
  double y_lim_physics[num_cuts][num_variables] = {{1,1,1,1,1,1,1,1,1,1,1,1},
						   {1,1,1,1,1,1,1,1,1,1,1,1},
						   {1,1,1,1,1,1,1,1,1,1,1,1},
						   {1,1,1,1,1,1,1,1,1,1,1,1},
						   {0.2,0.12,0.17,0.13,0.1,0.4,1,0.1,0.1,0.1,0.1,0.1}}; //y limit for area normalization plot 

  double y_lim_physics_xsec[num_cuts][num_variables] = {{1,1,1,1,1,1,1,1,1,1,1,1},
							{1,1,1,1,1,1,1,1,1,1,1,1},
							{1,1,1,1,1,1,1,1,1,1,1,1},
							{1,1,1,1,1,1,1,1,1,1,1,1},
							{0.2,0.12,0.17,0.13,0.1,0.4,1,0.1,0.1,0.1,0.1,0.1}}; //y limit for xsec plot

}; //end of class definition

//[area_normalized_plot]
// Plots histograms from the 4 samples. Samples have already be area normalized by GENIE_selection.C
// h_empirical, h_nieves, h_susa: Empirical, Nieves, and SuSA predictions for specific variable
// y_lim: y limit of the plot
// title-> title for plot
// path: where to save the png and pdf
// variable: name used to save the plot
void plotting::area_normalized_plot(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, double y_lim,const char* title, const char* path, const char* variable){

  TCanvas* canv_area = new TCanvas("canv_area_normalized","canv_area_normalize",2000,1500);
  //canv_area->cd();
  canv_area->SetRightMargin(0.09);
  canv_area->SetLeftMargin(0.15);
  canv_area->SetBottomMargin(0.15);
  
  h_empirical->Draw("hist"); //empiricial + lwellyn smith
  h_empirical->SetLineColor(kGreen+2);
  h_empirical->SetLineWidth(4);

  h_empirical->SetTitle(Form("%s %s","",title)); //titles_cuts[i]    
  h_empirical->GetXaxis()->SetTitle(Form("%s",title));
  h_empirical->GetXaxis()->SetTitleSize(0.035);
  h_empirical->GetYaxis()->SetTitle("Fractional Number of Events (%)"); //Counts
  h_empirical->GetYaxis()->SetTitleSize(0.035);
  h_empirical->SetMaximum(y_lim); //0.015//10000,//1900

  h_nieves->Draw("histSAME"); //nieves
  h_nieves->SetLineColor(kBlue);
  h_nieves->SetLineWidth(4);

  h_susa->Draw("histSAME"); //susav2
  h_susa->SetLineColor(kRed);
  h_susa->SetLineWidth(4);

  //h_GCF->Draw("histSAME"); //gcf
  //h_GCF->SetLineColor(kYellow+3);
  //h_GCF->SetLineWidth(4);

  TLatex *t_area = new TLatex(); //T Latex stuff                                                                                                                                                                                              
  t_area->SetNDC();
  t_area->SetTextAlign(22);
  t_area->DrawLatex(0.77,0.88,Form("%s",prelim));
  TLegend* legend_area = new TLegend(0.5,0.5,0.8,0.83,"");
  legend_area->AddEntry(h_empirical,"#splitline{Empirical MEC + Lwellyn Smith QE}{with FSI (GENIE hA2018 Model)}","L");
  legend_area->AddEntry(h_nieves,"#splitline{Nieves MEC & QE with FSI}{(GENIE hA2018 Model)}","L");
  legend_area->AddEntry(h_susa,"#splitline{SuSav2 MEC + QE with FSI}{(GENIE hA2018 Model)}","L");
  //legend_area->AddEntry(h_GCF,"#splitline{GCF with FSI}{(GENIE hA2018 Model)}","L");
  legend_area->SetBorderSize(0);
  legend_area->SetTextSize(0.03);
  legend_area->SetFillColor(0);
  legend_area->Draw("same");
  canv_area->Print(Form("%s%s_area_norm.png",path,variable));
  canv_area->Print(Form("%s%sarea_norm.pdf",path,variable));	

}

//Uses the method outlined in https://arxiv.org/pdf/2101.11867.pdf (Pg. 60, Eq 100) to maake the cross-section plots
//Method: (dsigma/dT)_j = (sigma * n_j)/(N* delta T_j)
//Method (in words): single differential cross-section as a function of T in bin j = (integrated cross-section * number of events in bin j)/(Number of simulated events * width of bin j)
//Allows for a direct comparison of the different generators on the CC1mu2p cross-section
//INPUTS:
//
void plotting::cross_section_plot(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF, double y_lim,const char* title, const char* path, const char* variable){

  //here are the sigmas. Taken from the GENIE splines                                                                                                                                                                                                                                                  
  double sigma_empirical = 3.02249 * 1E-37;
  double sigma_nieves = 27.1682 * 1E-38;
  double sigma_susa = 38.1417 * 1E-38;

  //N is the total number of events generated in each sample
  double N_empirical = 4000000;
  double N_nieves = 2400000;
  double N_susa = 3800000;
  //double N_GCF = 500000;

  //Here is were we calculate the dsigma/dx and the SD of dsigma/dx                                                                                                                                                                                                                                    
  double n_bins = h_empirical->GetNbinsX();

  for(int i=1; i < n_bins+1; i++){
    double delta_x = h_empirical->GetBinWidth(i);

    //empirical                                                                                                                                                                                                                                                                                        
    double n_empirical = h_empirical->GetBinContent(i) * N_empirical;
    double value_empirical = (sigma_empirical * n_empirical)/( N_empirical * delta_x);
    double SD_empirical = (sigma_empirical)/(delta_x*N_empirical) * std::sqrt(((N_empirical - n_empirical)*n_empirical)/(N_empirical));
    h_empirical->SetBinContent(i,value_empirical);
    h_empirical->SetBinError(i,SD_empirical);
    std::cout<<"EMPIRICAL value: "<<value_empirical<<std::endl;

    //nieves                                                                                                                                                                                                                                                                                          
    double n_nieves = h_nieves->GetBinContent(i);
    double value_nieves = (sigma_nieves * n_nieves)/(N_nieves * delta_x);
    double SD_nieves = (sigma_nieves)/(delta_x*N_nieves) * std::sqrt(((N_nieves - n_nieves)*n_nieves)/(N_nieves));
    h_nieves->SetBinContent(i,value_nieves);
    h_nieves->SetBinError(i,SD_nieves);
    std::cout<<"NIEVES value: "<<value_nieves<<std::endl;

    //susa                                                                                                                                                                                                                                                                                             
    double n_susa = h_susa->GetBinContent(i);
    double value_susa = (sigma_susa * n_susa)/(N_susa * delta_x);
    double SD_susa = (sigma_susa)/(delta_x*N_susa) * std::sqrt(((N_susa - n_susa)*n_susa)/(N_susa));
    h_susa->SetBinContent(i,value_susa);
    h_susa->SetBinError(i,SD_susa);
    std::cout<<"SUSA value: "<<value_susa<<std::endl;

    /*//Dealing with the GCF is a bit trickier                                                                                                                                                                                                                                                              //Each event in the GCF is assigned a weight which modifies the CCQE differential cross section                                                                                                                                                                                                    
    double n_GCF = h_GCF->GetBinContent(i);
    double value_GCF = (n_GCF)/(N_GCF*delta_x);
    double SD_GCF = (1)/(delta_x*N_GCF) * std::sqrt(((N_GCF - n_GCF)*n_GCF)/(N_GCF));
    h_GCF->SetBinContent(i,value_GCF);
    h_GCF->SetBinError(i,SD_GCF);
    std::cout<<"GCF value: "<<value_GCF<<std::endl;*/

  } //end of for loop 

  TCanvas* canv_xsec = new TCanvas("canv_xsec","canv_xsec",2000,1500);
  //canv_xsec->cd();
  canv_xsec->SetRightMargin(0.09);
  canv_xsec->SetLeftMargin(0.15);
  canv_xsec->SetBottomMargin(0.15);
  
  h_empirical->Scale(1/1E-38);
  h_empirical->Draw("hist"); //empiricial + lwellyn smith
  h_empirical->SetLineColor(kGreen+2);
  h_empirical->SetLineWidth(4);
  h_empirical->GetXaxis()->SetTitle(Form("%s",title));
  h_empirical->GetYaxis()->SetTitle("Differential Cross-Section [10^{-38} cm^{2} / Argon]");
  h_empirical->GetXaxis()->SetTitleSize(0.035);
  h_empirical->GetYaxis()->SetTitleSize(0.035);
  h_empirical->SetMaximum(y_lim); //0.015//10000,//1900
  h_empirical->SetTitle(Form("%s %s","",title)); //titles_cuts[i]

  h_nieves->Scale(1/1E-38);
  h_nieves->Draw("histSAME"); //nieves
  h_nieves->SetLineColor(kBlue);
  h_nieves->SetLineWidth(4);

  h_susa->Scale(1/1E-38);
  h_susa->Draw("histSAME"); //susav2
  h_susa->SetLineColor(kRed);
  h_susa->SetLineWidth(4);

  //h_GCF->Scale(1/1E-38);
  //h_GCF->Draw("histSAME"); //gcf
  //h_GCF->SetLineColor(kYellow+3);
  //h_GCF->SetLineWidth(4);

  TLatex *t_xsec = new TLatex(); //T Latex stuff                                                                                                                                                                                              
  t_xsec->SetNDC();
  t_xsec->SetTextAlign(22);
  t_xsec->DrawLatex(0.77,0.88,Form("%s",prelim));
  TLegend* legend_xsec = new TLegend(0.5,0.5,0.8,0.83,"");
  legend_xsec->AddEntry(h_empirical,"#splitline{Empirical MEC + Lwellyn Smith QE}{with FSI (GENIE hA2018 Model)}","L");
  legend_xsec->AddEntry(h_nieves,"#splitline{Nieves MEC & QE with FSI}{(GENIE hA2018 Model)}","L");
  legend_xsec->AddEntry(h_susa,"#splitline{SuSav2 MEC + QE with FSI}{(GENIE hA2018 Model)}","L");
  //legend_xsec->AddEntry(h_GCF,"#splitline{GCF with FSI}{(GENIE hA2018 Model)}","L");
  legend_xsec->SetBorderSize(0);
  legend_xsec->SetTextSize(0.03);
  legend_xsec->SetFillColor(0);
  legend_xsec->Draw("same");
  canv_xsec->Print(Form("%s%s_xsec.png",path,variable));
  canv_xsec->Print(Form("%s%s_xsec.pdf",path,variable));	

}//end of cross section plot

void plotting::truth_pred_plot(){
}//end of truth prediction plot

//Used to make the pmiss_plot() from Raquel's old technote. No longer needed
///////////////////////////////////////////////////////////////////////////
/*void plotting::pmiss_plot(){

    TCanvas* fig0 = new TCanvas("fig0","fig0",3400,800);
  fig0->Divide(3,1);
  fig0->cd(1);
  h_pmissT->Draw("e1");
  h_pmissT->SetLineColor(kBlue-7);
  h_pmissT->SetLineWidth(4);
  h_pmissT->SetTitle("MCC9 GENIE: Missing Transverse Momentum");
  h_pmissT->GetXaxis()->SetTitle("P_{miss}^{T} (GeV/c)");
  h_pmissT->GetYaxis()->SetTitle("Counts");
  TLine* a0x = new TLine(0.259,0,0.259,0.02);
  a0x->Draw("same");
  a0x->SetLineColor(kBlack);
  a0x->SetLineWidth(4);
  t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
  fig0->cd(2);
  h_Enu->Draw("E1");
  h_Enu->SetLineColor(kBlue-7);
  h_Enu->SetLineWidth(4);
  h_Enu->SetTitle("MCC9 GENIE:Energy of the Neutrino");
  h_Enu->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h_Enu->GetYaxis()->SetTitle("Counts");
  t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
  fig0->cd(3);
  h_pMiss->Draw("E1");
  h_pMiss->SetLineColor(kBlue-7);
  h_pMiss->SetLineWidth(4);
  h_pMiss->SetTitle("MCC9 GENIE:Struck Neutron Momentum");
  h_pMiss->GetXaxis()->SetTitle("Momentum of Struck Neutron (GeV/c)");
  h_pMiss->GetYaxis()->SetTitle("Counts");
  TLine* a0y = new TLine(0.259,0,0.259,0.027);
  a0y->Draw("same");
  a0y->SetLineColor(kBlack);
  a0y->SetLineWidth(4);
  t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
  fig0->cd();
  fig0->Print("images/raquel_24673.png");
  fig0->Print("images/raquel_24673.pdf");

  TCanvas* fig1 = new TCanvas("fig1","fig1",2000,1500);
  fig1->cd();
  h_2D_p1vp2->Draw("colz");
  h_2D_p1vp2->GetXaxis()->SetTitle("P_{Recoil} (GeV/c)");
  h_2D_p1vp2->GetYaxis()->SetTitle("P_{Leading} (GeV/c)");
  h_2D_p1vp2->SetTitle("");
  TLine* a1x = new TLine(0.259,0,0.259,3);
  a1x->Draw("same");
  a1x->SetLineColor(kBlue);
  a1x->SetLineWidth(4);
  TLine* a1y = new TLine(0,0.259,3,0.259);
  a1y->Draw("same");
  a1y->SetLineColor(kBlue);
  a1y->SetLineWidth(4);
  fig1->Print("images/raquel_Fig1.png");
  fig1->Print("images/raquel_Fig1.pdf");

  TCanvas* fig2 = new TCanvas("fig2","fig2",2000,1500);
  fig2->cd();
  h_cos_gamma_lab->Draw("E1");
  h_cos_gamma_lab->GetXaxis()->SetTitle("cos(#gamma_{Lab})");
  h_cos_gamma_lab->GetYaxis()->SetTitle("Counts");
  h_cos_gamma_lab->SetTitle("MCC9 GENIE: cos(#gamma_{Lab})");
  fig2->Print("images/raquel_Fig2.png");
  fig2->Print("images/raquel_Fig2.pdf");

  TCanvas* fig3 = new TCanvas("fig3","fig3",2000,1500);
  fig3->cd();
  h_2D_cos_gamma_cm_vs_p_missT->Draw("colz");
  h_2D_cos_gamma_cm_vs_p_missT->GetXaxis()->SetTitle("P_{miss}^{T} (GeV/c)");
  h_2D_cos_gamma_cm_vs_p_missT->GetYaxis()->SetTitle("cos(#gamma_{COM})");
  h_2D_cos_gamma_cm_vs_p_missT->SetTitle("");
  TLine* a3 = new TLine(0.3,-1,0.3,1);
  a3->Draw("same");
  a3->SetLineColor(kBlue);
  a3->SetLineWidth(4);
  fig3->Print("images/raquel_Fig3.png");
  fig3->Print("images/raquel_Fig3.pdf");

  TCanvas* fig4 = new TCanvas("fig4","fig4",2000,1500);
  fig4->cd();
  h_cos_gamma_cm_with_cut->Draw("e1");
  h_cos_gamma_cm_with_cut->GetXaxis()->SetTitle("cos(#gamma_{COM})");
  h_cos_gamma_cm_with_cut->GetYaxis()->SetTitle("Counts");
  h_cos_gamma_cm_with_cut->SetMaximum(12000);
  h_cos_gamma_cm_with_cut->SetTitle("MCC9 GENIE:cos(#gamma_{COM}) with Cut on P_{miss}^{T}");
  fig4->Print("images/raquel_Fig4.png");
  fig4->Print("images/raquel_Fig4.pdf");
} //end of pmiss_plot */
