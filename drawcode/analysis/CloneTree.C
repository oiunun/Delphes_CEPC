#ifdef __CLING__
R__LOAD_LIBRARY(../../libDelphes.so)
#endif
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
#include "TMath.h"
#include <string>
#include "TFile.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TClonesArray.h"

double Inv_Mass(Jet *jet,TLorentzVector &P_pi,TLorentzVector &P_em,TLorentzVector &P_ep);
void CloneTree()
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("../../../rootfile/mal/output_pythia6_PID/bb_*.root");
  TFile *f = new TFile("bgana.root","RECREATE");
  Int_t id_entry;
  TTree * t1 = new TTree("tree","id");
  t1->Branch("id_entry",&id_entry);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchPFC = treeReader->UseBranch("ParticleFlowCandidate");
  TH1F *hist_bg = new TH1F("bg", "bg", 100,4,6.6);

  double mass_b=5.27934;
  Long64_t numberOfEntries = treeReader->GetEntries();
  for(Long64_t entry = 0; entry <numberOfEntries ; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    if(branchJet->GetEntriesFast()!= 2) continue;
    Double_t Mass[2]={0.0,0.0};
    Double_t mass=0.0;
    TLorentzVector P_pi,P_em,P_ep;
    TLorentzVector P_pi_1,P_em_1,P_ep_1;
    TLorentzVector P_pi_2,P_em_2,P_ep_2;
    for(Int_t i=0;i<2;i++){
      Jet* jet=(Jet*) branchJet->At(i);
      if(i==0)  Mass[i]=Inv_Mass(jet,P_pi_1,P_em_1,P_ep_1);
      else  Mass[i]=Inv_Mass(jet,P_pi_2,P_em_2,P_ep_2);
    }
    if(abs(Mass[0]-mass_b)<abs(Mass[1]-mass_b)){
      mass=Mass[0];
      P_pi=P_pi_1;
      P_em=P_em_1;
      P_ep=P_ep_1;
    }
    else{
      mass=Mass[1];
      P_pi=P_pi_2;
      P_em=P_em_2;
      P_ep=P_ep_2;
    }
   if(P_pi.Pt()<1 || P_em.Pt()<1 || P_ep.Pt()<1) continue;
   if(mass<4 || mass>6.6) continue;
   hist_bg->Fill(mass);
   id_entry=entry;
   t1->Fill();
  } 
  t1->Write();
  f->Close();
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  gStyle->SetOptStat("im");
  hist_bg->Draw();
  c1->Print("fig/B_e1e1pi/hist_bg.png"); 
  c1->Print("fig/B_e1e1pi/hist_bg.pdf"); 
  c1->Clear();

  delete c1;
  delete treeReader;
  
}

double Inv_Mass(Jet *jet,TLorentzVector &P_pi,TLorentzVector &P_em,TLorentzVector &P_ep){
    int n_Jetp = jet->Constituents.GetEntriesFast();
    if(n_Jetp==0) return 0;
    std::vector < Int_t > id_pi;
    std::vector < Int_t > id_ep;
    std::vector < Int_t > id_em;
    ParticleFlowCandidate *particleflowcandidate=NULL;
    for(int i=0;i<n_Jetp;i++){
       particleflowcandidate = (ParticleFlowCandidate*) jet->Constituents.At(i);
      if((particleflowcandidate->Eem+particleflowcandidate->Ehad)/particleflowcandidate->P>0.4 && (particleflowcandidate->Eem/particleflowcandidate->P)<0.8   && abs(particleflowcandidate->PID)==211 /* &&  particleflowcandidate->P>3 */){
        id_pi.push_back(i);
      }
      if(particleflowcandidate->Eem/particleflowcandidate->P>0.8  /* &&  particleflowcandidate->P>3 */  && particleflowcandidate->Charge == 1 ){
        id_ep.push_back(i);
      }
      if(particleflowcandidate->Eem/particleflowcandidate->P>0.8 /* && particleflowcandidate->P>3 */ && particleflowcandidate->Charge == -1){
        id_em.push_back(i);
      }
      
    }

    float inv_mass = 0;
    if(!id_em.empty() && !id_ep.empty() && !id_pi.empty()) {       
      vector<int>::iterator it_pi;
      vector<int>::iterator it_ep;
      vector<int>::iterator it_em;
      TLorentzVector P_pi_i; 
      TLorentzVector P_ep_i;
      TLorentzVector P_em_i;
      for(it_pi=id_pi.begin();it_pi!=id_pi.end();it_pi++){
        ParticleFlowCandidate* canpi = (ParticleFlowCandidate*) jet->Constituents.At(*it_pi);
        P_pi_i=canpi->P4();
        for(it_ep=id_ep.begin();it_ep!=id_ep.end();it_ep++){
         ParticleFlowCandidate* canep = (ParticleFlowCandidate*) jet->Constituents.At(*it_ep);
          P_ep_i = canep->P4();
          for(it_em=id_em.begin();it_em!=id_em.end();it_em++){
           ParticleFlowCandidate* canem = (ParticleFlowCandidate*) jet->Constituents.At(*it_em);
            P_em_i = canem->P4();
            Double_t mass = (P_pi_i+P_ep_i+P_em_i).M();
            if(abs(mass-5.27934)<abs(inv_mass-5.27934)) {
              inv_mass = mass;
	      P_pi=P_pi_i;
	      P_em=P_em_i;
	      P_ep=P_ep_i;
            }
          }
        }
      }
    }
    id_em.clear();
    id_ep.clear();
    id_pi.clear();
    return inv_mass;
}
