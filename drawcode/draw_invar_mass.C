/*
root -l examples/draw_invar_mass.C'("delphes_output.root")'
*/

//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "../classes/DelphesClasses.h"
#include "../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
//#endif

//------------------------------------------------------------------------------

void draw_invar_mass(const char * inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  TLorentzVector cms(0,0,0,240);
  // Book histograms
  TH1 *histMass_muon = new TH1F("mass", "M_{inv}(#mu_{1}, #mu_{2})", 100, 75, 110.);
  TH1 *histMass_recoil = new TH1F("mass_rec", "M_{recoil}", 100, 120, 140.);
  TH1 *histMass_rr = new TH1F("mass_rr" , "M_{inv}(r_{1}, r_{2})", 100, 120., 130.);
  TH1 *histMass_jet = new TH1F("mass_jet" , "M_{inv}(jet_{1}, jet_{2})", 100, 100., 150.);
   
// Int_t total_nomatch = 0,total_excessmatch = 0;
  
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    if(branchJet->GetEntries() > 1  && branchMuon->GetEntries() > 1 /* && branchPhoton->GetEntries() > 1 */)
    {
      Jet *jet1 = (Jet*) branchJet->At(0);
      Jet *jet2 = (Jet*) branchJet->At(1);
      

      histMass_jet->Fill((jet1->P4()+jet2->P4()).M());

    }

    Muon *mu1, *mu2;

    // If event contains at least 2 muons
    if(branchMuon->GetEntries() > 1 && branchJet->GetEntries() > 1)
    {
      // Take first two muons
      mu1 = (Muon *) branchMuon->At(0);

      mu2 = (Muon *) branchMuon->At(1);

      // Plot their recoil mass
      histMass_recoil->Fill(((cms-mu1->P4()) -(mu2->P4())).M());
      histMass_muon->Fill((mu1->P4() + mu2->P4()).M());
    }

    if(branchPhoton->GetEntries() > 1)
    {
      Photon *r1 = (Photon*) branchPhoton->At(0);
      Photon *r2 = (Photon*) branchPhoton->At(1);

      histMass_rr->Fill((r1->P4()+r2->P4()).M());

    }
   }
 
    TCanvas *c1=new TCanvas("c1","c1",700,700);
    histMass_jet->Draw();
    histMass_jet->Fit("gaus");
    c1->Print("fig/histMass_jet.pdf"); 
    c1->Print("fig/histMass_jet.png"); 
    std::cout<<"numberOfEntries = "<<numberOfEntries<<endl;
    
    
    // histMass_rr->Draw();

    // histMass_rr->Fit("gaus");
    
/* a simple code to realize the MuonFilter module , but you must select the ParticleFlowCandidate branch
   ,and it need to be improved.

    Muon *mu[5] ;  Muon *mu_pos[5];  Muon *mu_neg[5];
    Muon *mu1 ,*mu2;
    // If event contains at least 2 electrons
    if(branchMuon->GetEntries() > 1)
    {
      Int_t n1 = 0,n2 = 0,nmatch =0 ;
      if(branchMuon->GetEntries() > 2){
         Int_t n_mu = branchMuon->GetEntries();
         //read all muons
         for(Int_t i=0 ;i<n_mu ;i++)       {mu[i] = (Muon *) branchMuon->At(i);}
         //part mu+ and mu-
         for(Int_t i=0 ;i<n_mu ;i++){
           if(mu[i]->Charge > 0) {mu_pos[n1] = mu[i] ; n1++;}
           else {mu_neg[n2] = mu[i] ; n2++;}
         }
         //choose correct mu+ and mu-
         for(Int_t i=0 ;i<n1 ;i++){
           for(Int_t j=0 ;j<n2 ;j++){
             Float_t p1 = (mu_pos[i]->P4()).P();
             Float_t p2 = (mu_neg[j]->P4()).P();
             Float_t s = abs((mu_pos[i]->P4()+mu_neg[i]->P4()).M()-91.5)
             if(p1>18. && p2>18. && s<15.){
               mu1 = mu_pos[i];
               mu2 = mu_neg[j]; nmatch++;
             }
           }
         }
         if(nmatch ==0) total_nomatch++;
         if(nmatch >1) total_excessmatch++;
      }
      else {  
      // Take correct two muons
      mu1 = (Muon *) branchMuon->At(0);
      mu2 = (Muon *) branchMuon->At(1);
      nmatch++;
      }
      // Plot their invariant mass
      if(nmatch !=0)  histMass->Fill(((cms-mu1->P4()) - (mu2->P4())).M());
    }
  }
  std::cout << "total entries  "<<numberOfEntries << endl;
  std::cout << "matching entries  "<<histMass->GetEntries() << endl; 
  std::cout << "no matching two muons  "<<total_nomatch << endl;
  std::cout << "can not find the correct two muons "<<total_excessmatch <<endl;
  // Show resulting histograms
  histMass->Draw();
  histMass->Fit("gaus");*/
}

