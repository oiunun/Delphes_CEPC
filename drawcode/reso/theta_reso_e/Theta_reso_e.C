#define theta_reso_e_cxx
#include "theta_reso_e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>

void theta_reso_e::Loop()
{
//   In a ROOT session, you can do:
//      root> .L phi_reso_pt.C
//      root> phi_reso_pt t
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

  TH1F *t_gaus1[30];       //each bin histogram
  for(Int_t i=0;i<30;i++){ char t_hname[30];
    sprintf(t_hname,"hg1%02d",i);
    t_gaus1[i]=new TH1F(t_hname,"",200,-0.05,0.05);}

  // Loop over all events
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      for ( int j =0;  j<Photon_size; j++){

      Float_t truth_Theta = acos(Photon_Truth_CosTheta[j]);
      Float_t Theta = acos(Photon_CosTheta[j]);
      Float_t var = truth_Theta-Theta;
       if(Photon_E[j]>1.) {
        for ( int jp=0;jp<30;jp++){ if(Photon_E[j]>3.*jp && Photon_E[j]<=3.*jp+3.)
          {t_gaus1[jp]->Fill(var); } }
     } } 
  }

  double t_error1[30];
  for(Int_t i=0;i<30;i++){
    t_error1[i]=t_gaus1[i]->GetStdDev();
    delete t_gaus1[i];
   }

  TH1F *gaus1[30];       //each bin histogram
  for(Int_t i=0;i<30;i++){ char hname[30];
    sprintf(hname,"hg1%02d",i);
    gaus1[i]=new TH1F(hname,"",200,0-5*t_error1[i],0+5*t_error1[i]);}
  TH1F *pth1[30];        //each phi bin histogram
  for(Int_t i=0;i<30;i++){ char pthname[30];
    sprintf(pthname,"hpt1%02d",i);
    pth1[i]=new TH1F(pthname,"",100,3.*i,3.*i+3.);}


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      for ( int j =0;  j<Photon_size; j++){
      Float_t truth_Theta = acos(Photon_Truth_CosTheta[j]);
      Float_t Theta = acos(Photon_CosTheta[j]);
      Float_t var = truth_Theta-Theta;
       if(Photon_E[j]>1.) {
        for ( int jp=0;jp<30;jp++){ if(Photon_E[j]>3.*jp && Photon_E[j]<=3.*jp+3.)
          {gaus1[jp]->Fill(var); pth1[jp]->Fill(Photon_E[j]);} }
     }}
   }
    double error1[30];
    double pt1[30];
    for(Int_t i=0;i<30;i++){
    pt1[i]=pth1[i]->GetMean();
    error1[i]=gaus1[i]->GetStdDev();
   }
  
  // Show resulting histograms
  TGraph *g;
  g = new TGraph(30,pt1,error1);


  int mak[] = {20, 21, 22, 23, 24, 25, 26};
  int line[] = {1, 1, 1, 1, 4, 4, 4};  
  int col[] = {kAzure+2, kOrange+1, kGreen+2, kMagenta+2, kOrange+1, kGreen+2, kMagenta+2};  
  g->SetMarkerStyle(mak[0]);
  g->SetMarkerColor(col[0]);
  g->SetMarkerSize(0.7);
  g->SetLineColor(col[0]);
  g->SetLineStyle(line[0]);
  g->SetLineWidth(1);
  g->GetYaxis()->SetTitle("#sigma#theta(rad)");
  g->GetXaxis()->SetTitle("E(GeV)");
  g->GetYaxis()->SetRangeUser(2e-3,5e-3);
  TCanvas *c1=new TCanvas("c1","c1",1300,700);

  g->Draw("ALP");

  c1->Print("../fig/Theta_reso_e.pdf");
  c1->Print("../fig/Theta_reso_e.png");
  TFile *f=new TFile("../fig/root/his_Theta_e.root","recreate");
  for(int i = 0 ;i<30;i++){
    gaus1[i]->Write();
  }
  f->Close();

}
void Theta_reso_e(){
  theta_reso_e t;
  t.Loop();
}