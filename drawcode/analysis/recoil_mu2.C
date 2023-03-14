#ifdef __CLING__
R__LOAD_LIBRARY(../../libDelphes.so)
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
#endif
void recoil_mu2()
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("../../rootfile/e2e2h_X_1_9.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");

  GenParticle *particle;
  ParticleFlowCandidate *particleflowcandidate,*particlemu1,*particlemu2;
  Long64_t numberOfEntries = treeReader->GetEntries();
  TH1F* his_mass = new TH1F("h_mass","mass",100,120,140);

  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    vector<int> mu1;
    vector<int> mu2;
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);

      if(particle->PID==13) mu1.push_back(i);
      if(particle->PID==-13) mu2.push_back(i);
      
      }
    Double_t inv_mass = 0;
    if(!mu1.empty() &&!mu2.empty()){
      for(vector<int>::iterator it_mu1=mu1.begin();it_mu1!=mu1.end();it_mu1++){
        for(vector<int>::iterator it_mu2=mu2.begin();it_mu2!=mu2.end();it_mu2++){
          TLorentzVector P_mu1;
          TLorentzVector P_mu2;
          TLorentzVector cms(0,0,0,240);
          particlemu1= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(*it_mu1);
          particlemu2= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(*it_mu2);
          P_mu1.SetPtEtaPhiM(particlemu1->PT,particlemu1->Eta,particlemu1->Phi,0.10565);
          P_mu2.SetPtEtaPhiM(particlemu2->PT,particlemu2->Eta,particlemu2->Phi,0.10565);
          Double_t mass = (cms-P_mu1-P_mu2).M();
          if(abs(mass-125)<abs(inv_mass-125)) inv_mass = mass;
        }
      }
    }
    mu1.clear();
    mu2.clear();
    if(inv_mass>120 && inv_mass<140) his_mass->Fill(inv_mass);
  }
   TCanvas *c1=new TCanvas("c1","c1",700,700);
   his_mass->Draw();
   c1->Print("inv_mass.png");


}