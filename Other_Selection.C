////////////////////
//7/19/2022: Other_Selection.C
//Author: Samantha Sword-Fehlberg
//This was an effort to analyze CC1mu2p events from GiBBU and NEUT. 
//I never finished developing this cod as we didn't want to include it in the XSec measurement
//Others should feel free to work on this code.
///////////////////////
#define Other_Selection_cxx

void Other_Selection::Loop(){

  auto start = high_resolution_clock::now();  
  
  //Open output root file and resulting table
  ///////////////////////////////////////
  TFile* outputfile;
  if(fsi == true){
    outputfile = new TFile(Form("root_files/histograms/fsi/hists_%s_fsi.root",sample), "RECREATE");
  }else{
    outputfile = new TFile(Form("root_files/histograms/no_fsi/hists_%s_no_fsi.root",sample), "RECREATE");
  }
  std::ofstream result_table; //csv file containing the number of events after each cut
  std::ofstream interaction_table; //csv file containing the breakdown of events in terms of interaction modes
  result_table.open(Form("tables/%s_remaining_events.csv",sample));
  interaction_table.open(Form("tables/%s_interaction_modes.csv",sample));

  //Define all the Histograms
  ///////////////////////////
  histogram_funcs hist; //histograms_funcs.h
  hist.Define_Histograms(); //define all the histograms you want

  //random constants and things we will need later
  TVector3 temp;
  TVector3 vAll_protons_mom;
  TVector3 vProtons_mom;
  int muon_mom_mag = 0;
  int proton_mom_mag = 0;
  int pionpm_requirement = 0;
  int pion0_requirement = 0;
  int event_count = 0;
  double weight;
  double factor = 100 * std::pow(0.1973,3);
  double conversion = 1e-33;
  srand( unsigned(time(NULL)));

  //counters
  //number of selected event/channel                                                                             
  int cc0 = 0;
  int nc0 = 0;
  int qel0 = 0;
  int res0 = 0;
  int mec0 = 0;
  int coh0 = 0;
  int dis0 = 0;

  //Begin Loop over all the Events:
  /////////////////////////////////////
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry % 10000 == 0){
      std::cout<<" EVENT "<<jentry<<std::endl;
    }
     
     //Before we do anything make sure to set the correct weight:
     //////////////////////////////////////////////////////////
     if(strcmp(sample,"GCF_CCQE") == 0){
       weight = (wght/factor) * conversion;
     } else{
       weight = wght;
     }

     //Okay. So first we need to make sure we grab the right events:
     /////////////////////////////////////////////////////////////////
     TVector3 vMuon_all(pxl,pyl,pzl);
     std::vector<int> protons_id;
     std::vector<double> protons_mom;
     std::vector<std::pair<double,int>> zipped;

     //This portion will identify the number of threshold protons, charged pions, and neutral pions
     //It will also save the proton id's and momenta's in order to later identify the leading and recoil protons
     //0) Loop through all the final state particles (nf)
     //1) Calcualte the particle momentum (particle_mom), the pdg (pdg), and the rescattering number (rescatter)
     //2) If the particle meets the threshold requirements you can add it to the appropriate counter
     //NOTE: There is a loop if fsi== false and an else loop. The fsi == false loop is an attempt to study the events before fsi effects. Every particle has an additional
     //requirement that the rescatter constant == 1. This same requirement is not needed for the else (fsi turned on) loop.
     //Finding the indices of the leading and recoil protons
     ////////////////////////////////////////////////////////////
     int mc_n_threshold_proton = 0;
     int mc_n_threshold_pion0 = 0;
     int mc_n_threshold_pionpm = 0;

     for(int i=0; i < nf; i++){
       TVector3 particle_mom(pxf[i],pyf[i],pzf[i]);
       int pdg = pdgf[i];
       int rescatter = resc[i];

       //If we want to look at the events without FSI
       //NEED TO ENSURE TO ADD RESCATTER == 1 REQUIREMENT
       ///////////////////////////////////////////
       if(fsi == false){

	 //protons
	 if(std::abs(pdg) == 2212){
	   if(particle_mom.Mag() > PROTON_MOM_CUT_LOW && particle_mom.Mag() < PROTON_MOM_CUT_HIGH && rescatter == 1){
	     protons_id.push_back(i);
	     mc_n_threshold_proton++;
	   }

	   //charged pions
	 }else if (std::abs(pdg) == 211){
	   if(particle_mom.Mag() > CHARGED_PI_MOM_CUT && rescatter == 1){
	     mc_n_threshold_pionpm++;
	   }

	   //neutral pions
	 } else if (pdg == 111 && rescatter == 1){
	   mc_n_threshold_pion0++;
	 }

	 //If we want to look at the events with FSI
	 ////////////////////////////////////////////
       } else{

	 //protons
	 if(std::abs(pdg) == 2212){
	   if(particle_mom.Mag() > PROTON_MOM_CUT_LOW && particle_mom.Mag() < PROTON_MOM_CUT_HIGH){
	     protons_id.push_back(i);
	     mc_n_threshold_proton++;
	   }

	   //charged pions
	 }else if (std::abs(pdg) == 211){
	   if(particle_mom.Mag() > CHARGED_PI_MOM_CUT){
	     mc_n_threshold_pionpm++;
	   }

	   //neutral pions
	 } else if (pdg == 111){
	   mc_n_threshold_pion0++;
	 }

       } //end of else loop over FSI events
     } //end of loop over the number final state particles

     //Now we can begin the selection
     //0) We require the lepton to fall within our defined phase-space limits
     //1) We require there to be exactly 2 protons above threshold.
     //2) We also require there to be exactly 0 neutral pions and 0 charged pions above the threshold
     ////////////////////////////////
     
     if(vMuon_all.Mag() < MUON_MOM_CUT_LOW || vMuon_all.Mag() > MUON_MOM_CUT_HIGH) continue; //the lepton is above threshold
     muon_mom_mag++;

     if(mc_n_threshold_proton != 2) continue; //there are exactly two protons above threshold
     proton_mom_mag++;

     if(mc_n_threshold_pionpm != 0) continue;
     pionpm_requirement++;

     if(mc_n_threshold_pion0 != 0) continue; //no charged pions above threshold and no neutral pions
     pion0_requirement++;

     //Now we have to determine which proton is the leading (one with the most momentum) 
     //and the recoil proton (one with the least momentum)
     //0) Make a zipped object with the momentum of a proton and the proton id
     //1) Order the zipped object based on the momentum
     //2) Unpack the values into leading and recoil ids
     /////////////////////////////////////////////////////////////

     for(int j=0; j < protons_id.size(); j++){      
       int value=protons_id[j];
       vProtons_mom.SetXYZ(pxf[value],pyf[value],pzf[value]);
       protons_mom.push_back(vProtons_mom.Mag());
       zipped.push_back(std::make_pair(protons_mom[j],protons_id[j]));
     }

     std::sort(zipped.begin(), zipped.end(),greater());
      
     for(int k = 0; k < protons_id.size(); k++){
       protons_id[k] = zipped[k].second;
       protons_mom[k] = zipped[k].first;
     }

     int leading_id = protons_id[0];
     int recoil_id = protons_id[1];

     if(_debug) std::cout<<"Leading Proton ID: "<<leading_id<<std::endl;
     if(_debug) std::cout<<"Recoil Proton ID: "<<recoil_id<<std::endl;
   
     //Now we are going to define some quantities we need to fill histograms
     //Remember these are also defined in the .h. Any changes you make here must also go there
     ////////////////////////////////////////////////////////////////////////////////////////
     TVector3 vBeam(0.,0.,Ev); // z-direction is defined along the neutrino direction
     TVector3 vMuon(pxl,pyl,pzl); //muon momentum
     TVector3 vLead(pxf[leading_id],pyf[leading_id],pzf[leading_id]); //leading proton momentum
     TVector3 vRec(pxf[recoil_id],pyf[recoil_id],pzf[recoil_id]); //recoil proton momentum
     TVector3 pn_vec(pxn,pyn,pzn);
     double nu_true = Ev;

     if(_debug) std::cout<<"Leading Proton Momentum: "<<vLead.Mag()<<std::endl;
     if(_debug) std::cout<<"Recoil Proton Momentum: "<<vRec.Mag()<<std::endl;
   

     //Now to fill histograms
     //The histograms get filled based on the momentum threshold cuts
     //No cuts = 0
     //Muon momentum cut = 1
     //A cut on the missing momentum = 2
     //^^^(this is depreicated now. it was an attempt to make my plots look like Raquels)
     //Recoil Proton momentum cut = 3
     //Leading Proton momentum cut = 4
     //NOTE: AS OF 12/20/2021 I HAVE INCLUDED THE THRESHOLD REQUIREMENTS IN THE BEGINNING OF THIS CODE
     //I .E. THE PLOTS BELOW SHOULD NOT LOOK DIFFERENT FROM EACH OTHER
     ///////////////////////////////////////////////////////////////////

     //Fill Histograms before any cuts:
     hist.Fill_Histograms(0,vBeam, vMuon, vLead, vRec, nu_true,pn_vec,weight);
     n_events[0]++;

     //Cut Requiring |pMiss| > 300 MeV/c
     //if(vmiss.Mag() < 0.3) continue; // Placing a cut requiring |pMiss| > 300 MeV/c
     hist.Fill_Histograms(1,vBeam, vMuon, vLead, vRec, nu_true, pn_vec,weight);
     n_events[1]++;
      
     //Cut on the Muon Momentum
     if(vMuon.Mag() < MUON_MOM_CUT_LOW || vMuon.Mag() > MUON_MOM_CUT_HIGH) continue;
     hist.Fill_Histograms(2,vBeam, vMuon, vLead, vRec, nu_true, pn_vec,weight);
     n_events[2]++;
      
     //Cut on the Leading Proton Momentum                                                                                                                                                                                                           
     if(vLead.Mag() < PROTON_MOM_CUT_LOW || vLead.Mag() > PROTON_MOM_CUT_HIGH) continue;
     hist.Fill_Histograms(3,vBeam, vMuon, vLead, vRec, nu_true, pn_vec,weight);
     n_events[3]++;

     //Cut on the Recoil Proton Momentum
     if(vRec.Mag() < PROTON_MOM_CUT_LOW || vRec.Mag() > PROTON_MOM_CUT_HIGH) continue; 
     hist.Fill_Histograms(4,vBeam, vMuon, vLead, vRec, nu_true, pn_vec,weight);
     n_events[4]++;

     //Just some checks for myself
     /////////////////////////////
     if(_debug){
       std::cout<<"Number of threshold protons: "<<mc_n_threshold_proton<<std::endl;
       std::cout<<"Number of threshold pion0: "<<mc_n_threshold_pion0<<std::endl;
       std::cout<<"Number of threshold pionpm: "<<mc_n_threshold_pionpm<<std::endl;
       
       TVector3 particle_mom;
       std::cout<<"[END OF THE CUTS] Number of Final State particles: "<<nf<<std::endl;
       for(int i=0; i < nf; i++){
	 std::cout<<"[END OF THE CUTS] PDF of the Particles: "<<pdgf[i]<<std::endl;
	 particle_mom.SetXYZ(pxf[i],pyf[i],pzf[i]);
	 std::cout<<"[END OF THE CUTS] momentum magnitude of the particle: "<<particle_mom.Mag()<<std::endl;
	 std::cout<<"[END OF THE CUTS] rescattering  value: "<<resc[i]<<std::endl;
       }
     }

     /////////////
     //Fill Counters
     ////////////
     if(cc == true){
       cc0++;
     }
     if(nc == true){
       nc0++;
     }
     if(qel == true){
       qel0++;
     }
     if(mec == true){
       mec0++;
     }
     if(res == true){
       std::cout<<"Number of final state particles: "<<nf<<std::endl;
       for(int i=0; i < nf; i++){
	 TVector3 particle_mom(pxf[i],pyf[i],pzf[i]);
	 int pdg = pdgf[i];
	 if(std::abs(pdg) == 211 || pdg == 111){
	   std::cout<<"FUCK ME!"<<std::endl;
	 }
	 std::cout<<Form("Particle %d PDG: ",i)<<pdg<<", Momentum: "<<particle_mom.Mag()<<std::endl;
       }
       res0++;
     }
     if(dis == true){
       dis0++;
     }
     if(coh == true){
       coh0++;
     }


     ////////////
     //DEPRICATED
     /////////////

     //Do all the stuf for Raquel's Analysis
     //hist.compute_raquel(vBeam, vMuon, vLead, vRec, pn_vec,weight);
     //Remove events with P_missT greater than 300 MeV
     //if(p_missT >= 0.3) continue; //this is broken!!!!!!! I NEED TO FIX THIS!
     //n_pmissT++;
     //h_cos_gamma_cm_with_cut->Fill(cos_gamma_cm,weight); //Figure 4 in Raquel's Note

     //Make sure to clear vectors when you are done
     protons_id.clear();
     protons_mom.clear();

   } //end loop over the entries

   result_table << "Total Number of Events "<<nentries<<"\n";
   result_table << "Number of Events Remaining after Muon Threshold Cut "<< muon_mom_mag<<"\n";
   result_table << "Number of Events Remaining after Proton Threshold Cut "<< proton_mom_mag<<"\n";
   result_table << "Number of Events Remaining after Pionpm Requirement "<< pionpm_requirement<<"\n";
   result_table << "Number of Events Remaining after Pion0 Requirement "<< pion0_requirement<<"\n";
   result_table << "Number of Events Remaining after pMiss Cut "<< n_events[1]<< "\n";
   result_table << "Number of Events Remaining after Muon Cut "<< n_events[2]<< "\n";
   result_table << "Number of Events Remaining after Lead Cut "<< n_events[3]<< "\n";
   result_table << "Number of Events Remaining after Recoil Cut "<< n_events[4]<< "\n";
   result_table << "Number of Events from Raquel's bit of code "<< n_pmissT<< "\n";
   result_table.close();

   interaction_table <<"Initial Number of Events That were Reconstructed"<< n_events[4] <<"\n";
   interaction_table <<"Number of CC Events "<< cc0 <<" \n";
   interaction_table <<"Number of NC Events "<< nc0 <<" \n";
   interaction_table <<"Number of QEL Events "<< qel0 <<" \n";
   interaction_table <<"Number of RES Events "<< res0 <<" \n";
   interaction_table <<"Number of MEC Events "<< mec0 <<" \n";
   interaction_table <<"Number of COH Events "<< coh0 <<" \n";
   interaction_table <<"Number of DIS Events "<< dis0 <<" \n";
   interaction_table.close();

   std::cout<<"---------[RESULTS OF CUTS]---------"<<std::endl;
   std::cout<<"Total Number of Events: "<<nentries<<std::endl;
   std::cout<<"Number of Events Remaining after Muon Threshold Cut: "<<muon_mom_mag<<" ("<<float(100.*(float(muon_mom_mag)/float(nentries)))<<"%)"<<" This is a difference of "<<nentries-muon_mom_mag<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after 2 Proton within Threshold Cut: "<<proton_mom_mag<<" ("<<float(100.*(float(proton_mom_mag)/float(nentries)))<<"%)"<<" This is a difference of "<<nentries-proton_mom_mag<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Pionpm Requirements: "<<pionpm_requirement<<" ("<<float(100.*(float(pionpm_requirement)/float(nentries)))<<"%)"<<" This is a difference of "<<nentries-pionpm_requirement<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Pion0 Requirements: "<<pion0_requirement<<" ("<<float(100.*(float(pion0_requirement)/float(nentries)))<<"%)"<<" This is a difference of "<<nentries-pion0_requirement<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after pMiss Cut: "<<n_events[1]<<" ("<<float(100.*(float(n_events[1])/float(nentries)))<<"%)"<<" This is a difference of "<<nentries-n_events[1]<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Muon Cut: "<<n_events[2]<<" ("<<float(100.*(float(n_events[2])/float(nentries)))<<"%)"<<" This is a difference of "<<n_events[1]-n_events[2]<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Lead Cut: "<<n_events[3]<<" ("<<float(100.*(float(n_events[3])/float(nentries)))<<"%)"<<" This is a difference of "<<n_events[2]-n_events[3]<<" Events."<<std::endl;
   std::cout<<"Number of Events Remaining after Recoil Cut: "<<n_events[4]<<" ("<<float(100.*(float(n_events[4])/float(nentries)))<<"%)"<<" This is a difference of "<<n_events[3]-n_events[4]<<" Events."<<std::endl;
   std::cout<<"Number of Events from Raquel's bit of Code: "<<n_pmissT<<std::endl;
   std::cout<<"------------------"<<std::endl;

   std::cout<<"---------[INTERACTION MODES]---------"<<std::endl;
   std::cout << "Initial Number of Events That were Reconstructed: "<<n_events[4]<<std::endl;
   std::cout << "Number of CC Events: "<<cc0<<" Fraction of the Total: "<<float(100.*(float(cc0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout << "Number of NC Events: "<<nc0<<" Fraction of the Total: "<<float(100.*(float(nc0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout << "Number of QEL Events: "<<qel0<<" Fraction of the Total: "<<float(100.*(float(qel0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout << "Number of RES Events: "<<res0<<" Fraction of the Total: "<<float(100.*(float(res0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout << "Number of MEC Events: "<<mec0<<" Fraction of the Total: "<<float(100.*(float(mec0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout << "Number of COH Events: "<<coh0<<" Fraction of the Total: "<<float(100.*(float(coh0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout << "Number of DIS Events: "<<dis0<<" Fraction of the Total: "<<float(100.*(float(dis0)/float(n_events[4])))<<"%"<<std::endl;
   std::cout<<"------------------"<<std::endl;

   //Write the histograms and 
   //close the output root file
   outputfile->cd();
   hist.Write_Histograms();
   outputfile->Close();

   //Program Clock
   ////////////////////////////////
   std::cout<<"----PROGRAM COMPLETED----"<<std::endl;
   auto stop = high_resolution_clock::now();
   auto duration = duration_cast<minutes>(stop - start); 
   std::cout<<"Program Run Time: "<<duration.count()<<std::endl;
   std::cout<<"-----------------------"<<std::endl;

}
