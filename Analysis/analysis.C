#include "helper_funcs.h"
#include <ctime>
#include <string>

void analysis(){

  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);

  //Load in the histogram file
  std::vector<TFile*> files; //vector of all the files
  TFile* f1=new TFile("/uboone/data/users/sfehlber/MC_Studies/GENIE/histograms/fsi/hists_empirical_lwellyn_fsi.root");
  TFile* f2=new TFile("/uboone/data/users/sfehlber/MC_Studies/GENIE/histograms/fsi/hists_nieves_fsi.root");
  TFile* f3=new TFile("/uboone/data/users/sfehlber/MC_Studies/GENIE/histograms/fsi/hists_susav2_fsi.root");
  TFile* f4=new TFile("/uboone/data/users/sfehlber/MC_Studies/GENIE/histograms/fsi/hists_GCF_CCQE_fsi.root"); //don't know how to compare to this yet. 

  files.push_back(f1);
  files.push_back(f2);
  files.push_back(f3);
  //files.push_back(f4);
  
  //stuff for date and time
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;
  
  const char* pathname = Form("/uboone/data/users/sfehlber/MC_Studies/images/%d%d%d/",Month,Day,Year);
  string path(pathname);
  int dir_exists = dirExists(pathname);
  if(dir_exists == 0){
    mkdir(pathname,0777);
    std::cout<<"New Directory Succesfully Created"<<std::endl;
  } else if (dir_exists < 0){
    std::cout<<"An Error has Occured. Please Check the Input Path Name."<<std::endl;
  }else if(dir_exists > 0){
    std::cout<<"Directory Already Exists. Continuing with Analysis"<<std::endl;
    }
  
  //Definitions for all the histograms we are going to make
  //////////////////////////////////////////////////////////
  
  //General stuff
  TLatex *t = new TLatex(); //T Latex stuff                                                                                                                                                        
  t->SetNDC();
  t->SetTextAlign(22);

  const int num_mc = 3; //normally 4, but don't know how to compare GCF with the MEC stuff
  const int num_cuts = 5;
  const char* cuts[num_cuts] = {"_b4_cuts","_pmis_cut","_muon_cut","_rec_cut","_lead_cut"};
  const char* titles_cuts[num_cuts] = {"Before Cuts","After P_{Miss} Cut","After Muon Mom. Cut","After Leading Mom. Cut","After Recoil Mom. Cut"};

  //static const int num_cuts = 7;
  //const char* cuts[num_cuts] = {"_b4_cuts","_pmiss_cut","_muon_cut","_lead_cut","_rec_cut","_pionpm_cut","_pion0_cut"};

  //Individual particle plots
  const int num_var = 4;
  const char* var[num_var] = {"_mom","_E","_theta","_phi"};
  const char* titles_var[num_var] = {"Momentum (GeV/c)","Energy (GeV)","cos(#theta)","#phi (Rad)"};
  const char* in_prog= "#scale[0.6]{MicroBooNE In-Progress}";
  double lead_proton_ylim[] = {0.1,0.3,0.6,0.05};
  double rec_proton_ylim[] = {0.25,0.5,0.3,0.05};
  double muon_ylim[] = {0.09,0.12,0.6,0.05};
  TH1D* h_muon[num_mc][num_cuts][num_var];
  TH1D* h_recoil[num_mc][num_cuts][num_var];
  TH1D* h_leading[num_mc][num_cuts][num_var];
  TCanvas* canv[num_cuts][num_var];
  TLegend* leg_leading[num_cuts][num_var];
  TLegend* leg_recoil[num_cuts][num_var];
  TLegend* leg_muon[num_cuts][num_var];

  
  //crap ton of physics plots
  const int num_variables = 12;
  const char* physics_variables[num_variables] = {"_opening_angle_protons_lab","_opening_angle_protons_com",
						  "_opening_angle_mu_leading","_opening_angle_mu_both",
						  "_delta_PT","_delta_alphaT","_delta_phiT","_pn",
						  "_tot_E","_tot_E_minus_beam","_PT_squared","_nu_E"};//_true_minus_calculated_E_neutrino"};
  
  const char* titles_physics[num_variables] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})",
					       "cos(#gamma_{#mu, p_{L}})", "cos(#gamma_{#mu, p_{L}+P_{R}})",
					       "#deltaP_{T} (GeV/c)","#delta#alpha_{T} (Deg.)","#delta#phi_{T} (Deg.)","p_{n} (GeV/c)",
					       "Total Muon Energy + Total Kinetic Energy of Protons (GeV/c)","Total Energy Minus Beam Energy (MeV/c)","(P^{T}_{miss})^{2} (GeV^{2}/c^{2})","Energy of the Neutrino (GeV/c)"};//,"E^{True}_{#nu} - E^{Calculated}_{#nu} (GeV/c)"};
  TH1D* h_physics[num_mc][num_cuts][num_variables];
  TCanvas* canv_physics[num_cuts][num_variables];
  TLegend* legend_physics[num_cuts][num_variables];
  double y_lim_physics[num_cuts][num_variables] = {{1,1,1,1,1,1,1,1,1,1,1,1},
						   {1,1,1,1,1,1,1,1,1,1,1,1},
						   {1,1,1,1,1,1,1,1,1,1,1,1},
						   {1,1,1,1,1,1,1,1,1,1,1,1},
						   {0.2,0.12,0.17,0.13,0.1,0.4,1,0.1,0.1,0.1,0.1,0.1}};

 
  //Grab the histograms
  //////////////////////
  TH1D* h_Enu = (TH1D*)f1->Get("h_Enu");
  TH1D* h_pmissT = (TH1D*)f1->Get("h_pmissT");
  TH1D* h_pMiss = (TH1D*)f1->Get("h_pMiss");
  TH1D* h_cos_gamma_lab = (TH1D*)f1->Get("h_cos_gamma_lab");
  TH1D* h_cos_gamma_cm_with_cut = (TH1D*)f1->Get("h_cos_gamma_cm_with_cut");
  TH2D* h_2D_p1vp2 = (TH2D*)f1->Get("h_2D_p1vp2");
  TH2D* h_2D_cos_gamma_cm_vs_p_missT = (TH2D*)f1->Get("h_2D_cos_gamma_cm_v_p_missT"); 
  
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

  std::cout<<"Finished Grabbing HIstograms"<<std::endl;
  
  //Plotting Time!
  ////////////////
  /*  TCanvas* fig0 = new TCanvas("fig0","fig0",3400,800);
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
  */

  for(int i = 0; i < num_cuts; i++){
    std::cout<<"Value of i: "<<i<<std::endl;
    for(int j =0; j < num_variables; j++){
      std::cout<<"Value of j: "<<j<<std::endl;

      canv_physics[i][j] = new TCanvas(Form("c%s%s",cuts[i],physics_variables[j]),Form("c%s%s",cuts[i],physics_variables[j]),2000,1500);
      canv_physics[i][j]->cd();
      canv_physics[i][j]->SetRightMargin(0.09);
      canv_physics[i][j]->SetLeftMargin(0.15);
      canv_physics[i][j]->SetBottomMargin(0.15);
      
      h_physics[0][i][j]->Draw("hist"); //empiricial + lwellyn smith
      h_physics[1][i][j]->Draw("histSAME"); //nieves
      h_physics[2][i][j]->Draw("histSAME"); //susav2
      //h_physics[3][i][j]->Draw("histSAME"); //gcf
      h_physics[0][i][j]->SetLineColor(kGreen+2);
      h_physics[0][i][j]->SetLineWidth(4);
      h_physics[1][i][j]->SetLineColor(kBlue);
      h_physics[1][i][j]->SetLineWidth(4);
      h_physics[2][i][j]->SetLineColor(kRed);
      h_physics[2][i][j]->SetLineWidth(4);
      //h_physics[3][i][j]->SetLineColor(kYellow+3);
      //h_physics[3][i][j]->SetLineWidth(4);
      h_physics[0][i][j]->GetXaxis()->SetTitle(Form("%s",titles_physics[j]));
      h_physics[0][i][j]->GetYaxis()->SetTitle("Fractional Number of Events (%)"); //Counts
      h_physics[0][i][j]->GetXaxis()->SetTitleSize(0.035);
      h_physics[0][i][j]->GetYaxis()->SetTitleSize(0.035);
      h_physics[0][i][j]->SetMaximum(y_lim_physics[i][j]); //0.015//10000,//1900
      h_physics[0][i][j]->SetTitle(Form("%s %s","",titles_physics[j])); //titles_cuts[i]

      t->DrawLatex(0.77,0.88,"#scale[0.6]{MicroBooNE In-Progress}");
      legend_physics[i][j] = new TLegend(0.5,0.5,0.8,0.83,"");
      legend_physics[i][j]->AddEntry(h_physics[0][i][j],"#splitline{Empirical MEC + Lwellyn Smith QE}{with FSI (GENIE hA2018 Model)}","L");
      legend_physics[i][j]->AddEntry(h_physics[1][i][j],"#splitline{Nieves MEC & QE with FSI}{(GENIE hA2018 Model)}","L");
      legend_physics[i][j]->AddEntry(h_physics[2][i][j],"#splitline{SuSav2 MEC + QE with FSI}{(GENIE hA2018 Model)}","L");
      //legend_physics[i][j]->AddEntry(h_physics[3][i][j],"#splitline{GCF with FSI}{(GENIE hA2018 Model)}","L");
      legend_physics[i][j]->SetBorderSize(0);
      legend_physics[i][j]->SetTextSize(0.03);
      legend_physics[i][j]->SetFillColor(0);
      legend_physics[i][j]->Draw("same");
      canv_physics[i][j]->Print(Form("%s%s%s.png",path.c_str(),cuts[i],physics_variables[j]));
      //canv_physics[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),cuts[i],physics_variables[j]));	
	
    } //end loop over number of variables
      

    //Particle specific plots of momentum, energy, phi, and theta
    ////////////////////////////////////////////////////////////
    for(int j = 0; j < num_var; j++){
      
      canv[i][j] = new TCanvas(Form("c%s%s",cuts[i],var[j]),Form("c%s%s",cuts[i],var[j]),4200,800);
      canv[i][j]->Divide(3,1);
      canv[i][j]->cd(1);
      canv[i][j]->SetRightMargin(0.09);
      canv[i][j]->SetLeftMargin(0.15);
      canv[i][j]->SetBottomMargin(0.15);      
      h_leading[0][i][j]->Draw("hist");
      h_leading[1][i][j]->Draw("histSAME");
      h_leading[2][i][j]->Draw("histSAME");
      //h_leading[3][i][j]->Draw("histSAME"); //susav2
      h_leading[0][i][j]->SetLineColor(kGreen+2);
      h_leading[0][i][j]->SetLineWidth(4);
      h_leading[1][i][j]->SetLineColor(kBlue);
      h_leading[1][i][j]->SetLineWidth(4);
      h_leading[2][i][j]->SetLineColor(kRed);
      h_leading[2][i][j]->SetLineWidth(4);
      //h_leading[3][i][j]->SetLineColor(kYellow+3);
      //h_leading[3][i][j]->SetLineWidth(4);
      h_leading[0][i][j]->GetXaxis()->SetTitle(Form("%s",titles_var[j]));
      h_leading[0][i][j]->GetYaxis()->SetTitle("Fractional Number of Events (%)"); //Counts                                                                                                                                                   
      h_leading[0][i][j]->GetXaxis()->SetTitleSize(0.04);
      h_leading[0][i][j]->GetYaxis()->SetTitleSize(0.04);
      h_leading[0][i][j]->SetMaximum(lead_proton_ylim[j]);
      h_leading[0][i][j]->SetTitle("");
      leg_leading[i][j] = new TLegend(0.5,0.5,0.8,0.83,"");
      leg_leading[i][j]->AddEntry(h_leading[0][i][j],"#splitline{Empirical MEC + Lwellyn Smith QE}{with FSI (GENIE hA2018 Model)}","L");
      leg_leading[i][j]->AddEntry(h_leading[1][i][j],"#splitline{Nieves MEC & QE with FSI}{(GENIE hA2018 Model)}","L");
      leg_leading[i][j]->AddEntry(h_leading[2][i][j],"#splitline{SuSav2 MEC + QE with FSI}{(GENIE hA2018 Model)}","L");
      //leg_leading[i][j]->AddEntry(h_leading[3][i][j],"#splitline{GCF with FSI}{(GENIE hA2018 Model)}","L");
      leg_leading[i][j]->SetBorderSize(0);
      leg_leading[i][j]->SetTextSize(0.03);
      leg_leading[i][j]->SetFillColor(0);
      leg_leading[i][j]->Draw("same");
      t->DrawLatex(0.5,0.93,Form("#scale[0.8]{%s: Leading Proton}",titles_var[j]));
      t->DrawLatex(0.78,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

      canv[i][j]->cd(2);
      h_recoil[0][i][j]->Draw("hist");
      h_recoil[1][i][j]->Draw("histSAME");
      h_recoil[2][i][j]->Draw("histSAME");
      //h_recoil[3][i][j]->Draw("histSAME"); //susav2
      h_recoil[0][i][j]->SetLineColor(kGreen+2);
      h_recoil[0][i][j]->SetLineWidth(4);
      h_recoil[1][i][j]->SetLineColor(kBlue);
      h_recoil[1][i][j]->SetLineWidth(4);
      h_recoil[2][i][j]->SetLineColor(kRed);
      h_recoil[2][i][j]->SetLineWidth(4);
      //h_recoil[3][i][j]->SetLineColor(kYellow+3);
      //h_recoil[3][i][j]->SetLineWidth(4);
      h_recoil[0][i][j]->GetXaxis()->SetTitle(Form("%s",titles_var[j]));
      h_recoil[0][i][j]->GetYaxis()->SetTitle("Fractional Number of Events (%)"); //Counts                                                                                                                                                    
      h_recoil[0][i][j]->GetXaxis()->SetTitleSize(0.04);
      h_recoil[0][i][j]->GetYaxis()->SetTitleSize(0.04);
      h_recoil[0][i][j]->SetMaximum(rec_proton_ylim[j]);
      h_recoil[0][i][j]->SetTitle("");
      leg_recoil[i][j] = new TLegend(0.5,0.5,0.8,0.83,"");
      leg_recoil[i][j]->AddEntry(h_recoil[0][i][j],"#splitline{Empirical MEC + Lwellyn Smith QE}{with FSI (GENIE hA2018 Model)}","L");
      leg_recoil[i][j]->AddEntry(h_recoil[1][i][j],"#splitline{Nieves MEC & QE with FSI}{(GENIE hA2018 Model)}","L");
      leg_recoil[i][j]->AddEntry(h_recoil[2][i][j],"#splitline{SuSav2 MEC + QE with FSI}{(GENIE hA2018 Model)}","L");
      //leg_recoil[i][j]->AddEntry(h_recoil[3][i][j],"#splitline{GCF with FSI}{(GENIE hA2018 Model)}","L");
      leg_recoil[i][j]->SetBorderSize(0);
      leg_recoil[i][j]->SetTextSize(0.03);
      leg_recoil[i][j]->SetFillColor(0);
      leg_recoil[i][j]->Draw("same");
      t->DrawLatex(0.5,0.93,Form("#scale[0.8]{%s: Recoil Proton}",titles_var[j]));
      t->DrawLatex(0.78,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

      canv[i][j]->cd(3);
      h_muon[0][i][j]->Draw("hist");
      h_muon[1][i][j]->Draw("histSAME");
      h_muon[2][i][j]->Draw("histSAME");
      //h_muon[3][i][j]->Draw("histSAME"); //susav2
      h_muon[0][i][j]->SetLineColor(kGreen+2);
      h_muon[0][i][j]->SetLineWidth(4);
      h_muon[1][i][j]->SetLineColor(kBlue);
      h_muon[1][i][j]->SetLineWidth(4);
      h_muon[2][i][j]->SetLineColor(kRed);
      h_muon[2][i][j]->SetLineWidth(4);
      //h_muon[3][i][j]->SetLineColor(kYellow+3);
      //h_muon[3][i][j]->SetLineWidth(4);
      h_muon[0][i][j]->GetXaxis()->SetTitle(Form("%s",titles_var[j]));
      h_muon[0][i][j]->GetYaxis()->SetTitle("Fractional Number of Events (%)"); //Counts
      h_muon[0][i][j]->GetXaxis()->SetTitleSize(0.04);
      h_muon[0][i][j]->GetYaxis()->SetTitleSize(0.04);
      h_muon[0][i][j]->SetTitle("");
      h_muon[0][i][j]->SetMaximum(muon_ylim[j]);
      leg_muon[i][j] = new TLegend(0.5,0.5,0.8,0.83,"");
      leg_muon[i][j]->AddEntry(h_muon[0][i][j],"#splitline{Empirical MEC + Lwellyn Smith QE}{with FSI (GENIE hA2018 Model)}","L");
      leg_muon[i][j]->AddEntry(h_muon[1][i][j],"#splitline{Nieves MEC & QE with FSI}{(GENIE hA2018 Model)}","L");
      leg_muon[i][j]->AddEntry(h_muon[2][i][j],"#splitline{SuSav2 MEC + QE with FSI}{(GENIE hA2018 Model)}","L");
      //leg_muon[i][j]->AddEntry(h_muon[3][i][j],"#splitline{GCF with FSI}{(GENIE hA2018 Model)}","L");
      leg_muon[i][j]->SetBorderSize(0);
      leg_muon[i][j]->SetTextSize(0.03);
      leg_muon[i][j]->SetFillColor(0);
      leg_muon[i][j]->Draw("same");
      t->DrawLatex(0.5,0.93,Form("#scale[0.8]{%s: Muon}",titles_var[j]));
      t->DrawLatex(0.78,0.88,"#scale[0.6]{MicroBooNE In-Progress}");

      canv[i][j]->cd();
      //t->DrawLatex(0.5,0.98,Form("#scale[0.95]{GCF Event Generator}"));      
      canv[i][j]->Print(Form("%s%s%s.png",path.c_str(),cuts[i],var[j]));
      //canv[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),cuts[i],var[j]));

    }
  }
  std::cout<<"----PROGRAM HAS FINISHED-----"<<std::endl;
}

