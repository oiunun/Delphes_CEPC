/*
root -l examples/draw_invar_mass.C'("delphes_output.root")'
*/
#ifdef __CLING__
R__LOAD_LIBRARY(../../libDelphes.so)
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
#endif

//------------------------------------------------------------------------------
void Eem_Ehad(ExRootTreeReader *treeReader);
void Truth_Mass(ExRootTreeReader *treeReader);

void B_e2e2pi()
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("../../rootfile/B_e2e2pi.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  // Book histograms

  // Eem_Ehad(treeReader);
  Truth_Mass(treeReader);
  delete treeReader;
  
}

void Eem_Ehad(ExRootTreeReader *treeReader){

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");

  GenParticle *particle;
  ParticleFlowCandidate *particleflowcandidate;
  Long64_t numberOfEntries = treeReader->GetEntries();
  TH1F *histEP_pi = new TH1F("EP_pi", "EP_pi", 100,0, 1.1);
  TH1F *histprobpi_mu = new TH1F("probpi_mu", "probpi_mu", 100,0, 1);
  TH1F *histEP_e = new TH1F("EP_e", "EP_e", 100,0, 2);
  TH1F *histEP_mu = new TH1F("EP_mu", "EP_mu", 100,0, 2);
  // TH1F *histprob_pi = new TH1F("prob_pi", "prob_pi", 100,0.998, 1);
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);

      if(abs(particle->PID)==211){
        Float_t EP_pi =(particleflowcandidate->Eem+particleflowcandidate->Ehad)/particleflowcandidate->P;
        if(EP_pi<1.1){
          histEP_pi->Fill(EP_pi);
        }
        // histprob_pi->Fill(particleflowcandidate->Prob_Pi);
      }
      if(abs(particle->PID)==13){
        histprobpi_mu->Fill(particleflowcandidate->Prob_Pi);
        Float_t EP_mu =(particleflowcandidate->Eem+particleflowcandidate->Ehad)/particleflowcandidate->P;
        histEP_mu->Fill(EP_mu);
      }
    }
  }
    TCanvas *c1=new TCanvas("c1","c1",900,700);
    histEP_mu->Draw();
    c1->Print("histEP_mu.png"); 
    // histprobpi_mu->Draw();
    // c1->Print("histprobpi_mu.png");
    delete c1;
}
void Truth_Mass(ExRootTreeReader *treeReader){
  
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");

  GenParticle *particle;
  GenParticle *par;
  ParticleFlowCandidate *particleflowcandidate;
  ParticleFlowCandidate *can1;
  ParticleFlowCandidate *can2;
  ParticleFlowCandidate *can3;
  GenParticle *p1;
  GenParticle *p2;
  GenParticle *p3;
  Long64_t numberOfEntries = treeReader->GetEntries();
  TH1F *hist_truth_mass = new TH1F("B_truth_mass", "B_truth_mass", 100,0,5.6);
  TH1F *hist_PID = new TH1F("PID", "PID", 100,-300,300);
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    Int_t id_B =0;
    Int_t id_pi =0;
    Int_t id_mu1 =0;
    Int_t id_mu2 =0;
    TLorentzVector P_pi(0,0,0,0); 
    TLorentzVector P_mu1(0,0,0,0);
    TLorentzVector P_mu2(0,0,0,0);
    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      par=(GenParticle* )branchParticle->At(i);
      if(abs(par->PID)==521){
        id_B=i;
      }
    }
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);

      if(particle->M1==id_B /* || particle->M1==id_B2 */ /* && particle->PT>0.5 */){
        switch(particle->PID){
          case 211  : 
            id_pi=i;
            P_pi=particleflowcandidate->P4();
            break;
          case -211  : 
            id_pi=i;
            P_pi=particleflowcandidate->P4();
            break;
          case 13 :
            id_mu1=i;
            P_mu1=particleflowcandidate->P4();
            break;
          case -13 :
            id_mu2=i;
            P_mu2=particleflowcandidate->P4();
            break;
        } 
      }
    }

    if(id_pi!=0 && id_mu1!=0 && id_mu2!=0){
        Float_t mass = (P_pi+P_mu1+P_mu2).M();
        // if(mass>5 && mass<5.6){
          hist_truth_mass->Fill(mass);
        // }
      }
    
  }
    TCanvas *c1=new TCanvas("c1","c1",700,700);
    // hist_truth_mass->Fit("gaus");
    hist_truth_mass->Draw();
    c1->Print("fig/B_e2e2pi/hist_truth_mass.png"); 
    // hist_PID->Draw();
    // c1->Print("PID.png");
    delete c1;
}