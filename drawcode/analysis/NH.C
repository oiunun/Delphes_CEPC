#ifdef __CLING__
R__LOAD_LIBRARY(../../libDelphes.so)
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
#include "TMath.h"
#endif


void NH()
{
    gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");

  // chain.Add("/cefs/higgs/wangshudong/delphes_run/output/E240.Pe2e2h_e2e2.e0.p0.whizard195/e2e2h_e2e2.e0.p0.00001.root");
  chain.Add("../../rootfile/e2e2e2e2.root");
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  // Book histograms

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");


  GenParticle *particle;
  ParticleFlowCandidate *particleflowcandidate;


  Long64_t numberOfEntries = treeReader->GetEntries();
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);

    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){
      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      // cout<<particleflowcandidate->truth_PID<<"    ";
      if(particleflowcandidate->truth_PID == 130){
        cout<<"particles number  ="<<particleflowcandidate->Particles.GetEntriesFast()<<endl;
        cout<<"particles pid "; 
        for(Long64_t j =0;j<particleflowcandidate->Particles.GetEntriesFast();j++){
          particle = (GenParticle*) particleflowcandidate->Particles.At(j);
          cout<<"  "<<particle->PID;
        }
        cout<<endl;
      cout<<"next NeutralHadron" <<endl;
      }
    }
    //     cout<<endl;
    //     cout<<"next event"<<endl;
    // for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
    //   // particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
    //   particle= (GenParticle*) branchParticle->At(i);
    //   if(particle->Status == 1 ) cout<<"  "<<particle->PID;
    // }
    cout<<endl;
    cout<<"next event"<<endl;;


  }
  delete treeReader;
}