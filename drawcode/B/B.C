#define b_cxx
#include "b.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>


void SetVector(TLorentzVector &V_P,Double_t Mass,Double_t Theta,Double_t Phi,Double_t P);


void b::Loop()
{
//   In a ROOT session, you can do:
//      root> .L b.C
//      root> b t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   TH1F* his_P_truth_pi = new TH1F("P_truth_pi","P_truth_pi",100,0,50);
   TH1F* his_P_mu = new TH1F("h_mu","h_mu",100,0,50);
   TH1F* his_P_pi = new TH1F("h_pi","h_pi",100,0,50);
   TH1F* his_mass = new TH1F("h_mass","mass",100,5.2,5.35);
   TH1F* his_truth_mass = new TH1F("h_truth_mass","truth_mass",100,5.2,5.35);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int id_B=0;
      // for(int i=0;i<Particle_size;i++){
      //   if(abs(Particle_PID[i])==521) id_B = i;
      //   if(id_B!=0 && Particle_M1[i]==id_B){
      //     if(abs(Particle_PID[i])==211) his_P_truth_pi->Fill(Particle_P[i]);
          
      //     if(abs(Particle_PID[i])==13) his_P_mu->Fill(Particle_P[i]);
      //   }
      // }
      
      vector<int> pi;
      vector<int> mu1;  //mu+
      vector<int> mu2;  //mu-

      for(int i=0;i<Track_size;i++){
        // if(abs(Track_Truth_PID[i])==211) his_P_pi->Fill(Track_Truth_P[i]);
        // if(abs(Track_Truth_PID[i])==13) his_P_mu->Fill(Track_Truth_P[i]);
        if(abs(Track_Truth_PID[i])==211 /* && Track_P[i]>5 */) pi.push_back(i);
        // if(Track_Truth_PID[i]==13) mu1.push_back(i);
        // if(Track_Truth_PID[i]==-13) mu2.push_back(i);
      }
      for(int i=0;i<Muon_size;i++){
        if(Muon_Charge[i]>0) mu1.push_back(i);
        if(Muon_Charge[i]<0) mu2.push_back(i);
      }
      // Double_t inv_mass = 0;
      // if(!pi.empty() && !mu1.empty() &&!mu2.empty()){
      //   vector<int>::iterator it_pi;
      //   vector<int>::iterator it_mu1;
      //   vector<int>::iterator it_mu2;
      //   for(it_pi=pi.begin();it_pi!=pi.end();it_pi++){
      //     for(it_mu1=mu1.begin();it_mu1!=mu1.end();it_mu1++){
      //       for(it_mu2=mu2.begin();it_mu2!=mu2.end();it_mu2++){
      //         TLorentzVector P_pi; 
      //         TLorentzVector P_mu1;
      //         TLorentzVector P_mu2;
      //         SetVector(P_pi,Track_Mass[*it_pi],acos(Track_Truth_CosTheta[*it_pi]),Track_Truth_Phi[*it_pi],Track_Truth_P[*it_pi]);
      //         SetVector(P_mu1,Track_Mass[*it_mu1],acos(Track_Truth_CosTheta[*it_mu1]),Track_Truth_Phi[*it_mu1],Track_Truth_P[*it_mu1]);
      //         SetVector(P_mu2,Track_Mass[*it_mu2],acos(Track_Truth_CosTheta[*it_mu2]),Track_Truth_Phi[*it_mu2],Track_Truth_P[*it_mu2]);
      //         Double_t mass = (P_pi+P_mu1+P_mu2).M();
      //         if(abs(mass-5.27934)<abs(inv_mass-5.27934)) inv_mass = mass;
      //       }
      //     }
      //   }  
      // }


      Double_t inv_mass = 0;
      if(!pi.empty() && !mu1.empty() &&!mu2.empty()){
        vector<int>::iterator it_pi;
        vector<int>::iterator it_mu1;
        vector<int>::iterator it_mu2;
        for(it_pi=pi.begin();it_pi!=pi.end();it_pi++){
          for(it_mu1=mu1.begin();it_mu1!=mu1.end();it_mu1++){
            for(it_mu2=mu2.begin();it_mu2!=mu2.end();it_mu2++){
              TLorentzVector P_pi; 
              TLorentzVector P_mu1;
              TLorentzVector P_mu2;
              SetVector(P_pi,Track_Mass[*it_pi],acos(Track_CosTheta[*it_pi]),Track_Phi[*it_pi],Track_P[*it_pi]);
              P_mu1.SetPtEtaPhiM(Muon_PT[*it_mu1],Muon_Eta[*it_mu1],Muon_Phi[*it_mu1],0.10565);
              P_mu2.SetPtEtaPhiM(Muon_PT[*it_mu2],Muon_Eta[*it_mu2],Muon_Phi[*it_mu2],0.10565);
              Double_t mass = (P_pi+P_mu1+P_mu2).M();
              if(abs(mass-5.27934)<abs(inv_mass-5.27934)) inv_mass = mass;
            }
          }
        }  
      }
      // if(inv_mass>5.2 && inv_mass<5.35 ){
      //   his_mass->Fill(inv_mass);
      // }
      if(inv_mass>5.2 && inv_mass<5.35){
        his_truth_mass->Fill(inv_mass);
      }
 
      // if (Cut(ientry) < 0) continue;
   }
   TCanvas *c1=new TCanvas("c1","c1",700,700);
  //  his_mass->Draw();
  //  c1->Print("inv_mass.png");
   his_truth_mass->Draw();
   c1->Print("truth_mass.png");
  //  his_P_pi->Draw();
  //  c1->Print("P_pi.png");
  //  his_P_truth_pi->Draw();
  //  c1->Print("truth_pi.png");
}
void SetVector(TLorentzVector &V_P,Double_t Mass,Double_t Theta,Double_t Phi,Double_t P){
  Double_t E = TMath::Sqrt(pow(Mass,2)+pow(P,2));
  Double_t Px = P*sin(Theta)*cos(Phi);
  Double_t Py = P*sin(Theta)*sin(Phi);
  Double_t Pz = P*cos(Theta);
  V_P.SetPxPyPzE(Px,Py,Pz,E);
}

void B(){
  b t;
  t.Loop();
  return 0;
}
