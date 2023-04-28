/*
root -l examples/draw_invar_mass.C'("delphes_output.root")'
*/
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

//------------------------------------------------------------------------------
void Eem_Ehad(ExRootTreeReader *treeReader);
void Signal_Mass(ExRootTreeReader *treeReader);
void Bg_ana(ExRootTreeReader *treeReadersig,ExRootTreeReader *treeReaderbg);
void Truth_Mass(ExRootTreeReader *treeReader);
void TrackEfficiency(ExRootTreeReader *treeReader);
void Performance_PID(ExRootTreeReader *treeReader);
void testcos(ExRootTreeReader *treeReader);
void Momentum(ExRootTreeReader *treeReader,ExRootTreeReader *treeReaderbg);
double Inv_Mass(Jet *jet,TLorentzVector &P_pi,TLorentzVector &P_em,TLorentzVector &P_ep,Double_t &l);

void B_e1e1pi()
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  TChain chainbg("Delphes");
  chain.Add("../../../rootfile/signal.root");
  chainbg.Add("../../../rootfile/mal/output_pythia6_PID/bb_0001.root");
  // chain.Add("../../rootfile/test_bb.root");
  // chain.Add("../../rootfile/qqh_X_1_9.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  ExRootTreeReader *treeReaderbg = new ExRootTreeReader(&chainbg);


  // Eem_Ehad(treeReader);
 // Signal_Mass(treeReaderbg);
  // Truth_Mass(treeReader);
   testcos(treeReader);
  // TrackEfficiency(treeReaderbg);
  // Performance_PID(treeReader);
 //  Momentum(treeReader,treeReaderbg);
 //  Bg_ana(treeReader,treeReaderbg);

  delete treeReader;
  delete treeReaderbg;
  
}

void Eem_Ehad(ExRootTreeReader *treeReader){

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");


  GenParticle *particle;
  ParticleFlowCandidate *particleflowcandidate;
  Track *eflowtrack;
  Long64_t numberOfEntries = treeReader->GetEntries();

  TH1F *histEem_mu = new TH1F("Eem_mu", "Eem/P for #mu", 100,0, 2);
  TH1F *histEhad_mu = new TH1F("Ehad_mu", "Ehad/P for #mu", 100,0, 2);
  TH1F *histEem_had = new TH1F("Eem_had", "Eem/P for hadron", 100,0, 2);
  TH1F *histEhad_had = new TH1F("Ehad_had", "Ehad/P for hadron", 100,0, 2);
  TH1F *histEem_e = new TH1F("Eem_e", "", 100,0, 2);
  TH1F *histEhad_e = new TH1F("Ehad_e", "", 100,0, 2);
  TH1F *histEem_pi = new TH1F("Eem_pi", "Eem/P for pi", 100,0, 2);
  TH1F *histEhad_pi = new TH1F("Ehad_pi", "Ehad/P for pi", 100,0, 2);
  TH1F *histEem_k = new TH1F("Eem_k", "Eem/P for k", 100,0, 2);
  TH1F *histEhad_k = new TH1F("Ehad_k", "Ehad/P for k", 100,0, 2);
  TH1F *histEem_p = new TH1F("Eem_p", "Eem/P for p", 100,0, 2);
  TH1F *histEhad_p = new TH1F("Ehad_p", "Ehad/P for p", 100,0, 2);

  TH1F *histEdepo_e = new TH1F("Edepo_e", "", 100,0, 2);
  TH1F *histEdepo_mu = new TH1F("Edepo_mu", "(Eem+Ehad)/P for #mu", 100,0, 2);
  TH1F *histEdepo_had = new TH1F("Edepo_had", "(Eem+Ehad)/P for hadron", 100,0, 2);
  TH1F *histEdepo_pi = new TH1F("Edepo_pi", "(Eem+Ehad)/P for pi", 100,0, 2);
  TH1F *histEdepo_k = new TH1F("Edepo_k", "(Eem+Ehad)/P for k", 100,0, 2);
  TH1F *histEdepo_p = new TH1F("Edepo_p", "(Eem+Ehad)/P for p", 100,0, 2);


  Double_t xEem[20],xEhad[20],xE_em[20],xE_had[20];
  TH1F *hisEem_reso[20];
  for(Int_t i=0;i<20;i++){ char hname[20];
  sprintf(hname,"hEem%02d",i);
  hisEem_reso[i]=new TH1F(hname,"",200,i*0.5,i*2);}
  TH1F *hisEhad_reso[20];
  for(Int_t i=0;i<20;i++){ char hname[20];
  sprintf(hname,"hEhad%02d",i);
  hisEhad_reso[i]=new TH1F(hname,"",200,i*0.5,i*2);}

  TH1F *hisE_em[20];        
  for(Int_t i=0;i<20;i++){ char Ehname[20];
    sprintf(Ehname,"hE_em%02d",i);
    hisE_em[i]=new TH1F(Ehname,"",100,i,i+1);}
  TH1F *hisE_had[20];        
  for(Int_t i=0;i<20;i++){ char Ehname[20];
    sprintf(Ehname,"hE_had%02d",i);
    hisE_had[i]=new TH1F(Ehname,"",100,i,i+1);}

  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);

      if(particleflowcandidate->truth_PID == 130) continue;
      if(abs(particle->PID)==13 ){
        
        Float_t Eem_mu =particleflowcandidate->Eem/particleflowcandidate->P;
        Float_t Ehad_mu =particleflowcandidate->Ehad/particleflowcandidate->P;
        histEem_mu->Fill(Eem_mu);
        histEhad_mu->Fill(Ehad_mu);
        histEdepo_mu->Fill(Ehad_mu+Eem_mu);
        for(int j=0;j<20;j++){
          if(particleflowcandidate->E*0.1<j+1 && particleflowcandidate->E*0.1>j){
            hisEem_reso[j]->Fill(particleflowcandidate->Eem);
            hisE_em[j]->Fill(particleflowcandidate->E*0.1);
          }
          if(particleflowcandidate->E*0.2<j+1 && particleflowcandidate->E*0.2>j){
            hisEhad_reso[j]->Fill(particleflowcandidate->Ehad);
            hisE_had[j]->Fill(particleflowcandidate->E*0.2);
          }
        }
      }
      if(abs(particle->PID)==11 ){
        Float_t Eem_e =particleflowcandidate->Eem/particleflowcandidate->P;
        Float_t Ehad_e =particleflowcandidate->Ehad/particleflowcandidate->P;
        histEem_e->Fill(Eem_e);
        histEhad_e->Fill(Ehad_e);
        histEdepo_e->Fill(Ehad_e+Eem_e);
        for(int j=0;j<20;j++){
          if(particleflowcandidate->E<j+1 && particleflowcandidate->E>j){
            hisEem_reso[j]->Fill(particleflowcandidate->Eem);
            hisE_em[j]->Fill(particleflowcandidate->E);
          }
        }
      }
      if(abs(particle->PID)==321 ){
        
        Float_t Eem_k =particleflowcandidate->Eem/particleflowcandidate->P;
        Float_t Ehad_k =particleflowcandidate->Ehad/particleflowcandidate->P;
        histEem_k->Fill(Eem_k);
        histEem_had->Fill(Eem_k);
        histEhad_k->Fill(Ehad_k);
        histEhad_had->Fill(Ehad_k);
        histEdepo_k->Fill(Ehad_k+Eem_k);
        histEdepo_had->Fill(Ehad_k+Eem_k);
        for(int j=0;j<20;j++){
          if(particleflowcandidate->E*0.3<j+1 && particleflowcandidate->E*0.3>j){
            hisEem_reso[j]->Fill(particleflowcandidate->Eem);
            hisE_em[j]->Fill(particleflowcandidate->E*0.3);
          }
          if(particleflowcandidate->E*0.7<j+1 && particleflowcandidate->E*0.7>j){
            hisEhad_reso[j]->Fill(particleflowcandidate->Ehad);
            hisE_had[j]->Fill(particleflowcandidate->E*0.7);
          }
        }
      }
      if(abs(particle->PID)==211 ){
        
        Float_t Eem_pi =particleflowcandidate->Eem/particleflowcandidate->P;
        Float_t Ehad_pi =particleflowcandidate->Ehad/particleflowcandidate->P;
        histEem_pi->Fill(Eem_pi);
        histEem_had->Fill(Eem_pi);
        histEhad_pi->Fill(Ehad_pi);
        histEhad_had->Fill(Ehad_pi);
        histEdepo_pi->Fill(Ehad_pi+Eem_pi);
        histEdepo_had->Fill(Ehad_pi+Eem_pi);
        for(int j=0;j<20;j++){
          if(particleflowcandidate->E*0.3<j+1 && particleflowcandidate->E*0.3>j){
            hisEem_reso[j]->Fill(particleflowcandidate->Eem);
            hisE_em[j]->Fill(particleflowcandidate->E*0.3);
          }
          if(particleflowcandidate->E*0.7<j+1 && particleflowcandidate->E*0.7>j){
            hisEhad_reso[j]->Fill(particleflowcandidate->Ehad);
            hisE_had[j]->Fill(particleflowcandidate->E*0.7);
          }
        }
      }
      if(abs(particle->PID)==2212 ){
        
        Float_t Eem_p =particleflowcandidate->Eem/particleflowcandidate->P;
        Float_t Ehad_p =particleflowcandidate->Ehad/particleflowcandidate->P;
        histEem_p->Fill(Eem_p);
        histEem_had->Fill(Eem_p);
        histEhad_p->Fill(Ehad_p);
        histEhad_had->Fill(Ehad_p);
        histEdepo_p->Fill(Ehad_p+Eem_p);
        histEdepo_had->Fill(Ehad_p+Eem_p);
        for(int j=0;j<20;j++){
          if(particleflowcandidate->E*0.3<j+1 && particleflowcandidate->E*0.3>j){
            hisEem_reso[j]->Fill(particleflowcandidate->Eem);
            hisE_em[j]->Fill(particleflowcandidate->E*0.3);
          }
          if(particleflowcandidate->E*0.7<j+1 && particleflowcandidate->E*0.7>j){
            hisEhad_reso[j]->Fill(particleflowcandidate->Ehad);
            hisE_had[j]->Fill(particleflowcandidate->E*0.7);
          }
        }
      }
    }
  }
  TFile *f=new TFile("fig/B_e1e1pi/Edepo_reso.root","recreate");
  for(int i=0;i<20;i++){
    hisEem_reso[i]->Write();
    hisEhad_reso[i]->Write();
    xE_em[i] = hisE_em[i]->GetMean();
    xE_had[i] = hisE_had[i]->GetMean();
    xEem[i] = hisEem_reso[i]->GetStdDev()/xE_em[i];
    xEhad[i] = hisEhad_reso[i]->GetStdDev()/xE_had[i];

  }

  TCanvas *c1=new TCanvas("c1","c1",1500,550);
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(kFALSE);
  c1->Divide(3,1);
  c1->cd(1);
  histEem_e->SetLineColor(kRed);
  histEem_mu->SetLineColor(kViolet);
  // histEem_pi->SetLineColor(kGreen);
  // histEem_k->SetLineColor(kViolet );
  
  // histEem_mu->Scale(1.0/histEem_mu->Integral());
  // histEem_e->Scale(1.0/histEem_e->Integral());
  // histEem_had->Scale(1.0/histEem_had->Integral());

  histEem_e->GetYaxis()->SetRangeUser(0,17000);
  // histEem_mu->GetYaxis()->SetRangeUser(0,1);
  // histEem_had->GetYaxis()->SetRangeUser(0,1);
  // histEem_pi->GetYaxis()->SetRangeUser(0,10000);
  // histEem_k->GetYaxis()->SetRangeUser(0,40000);
  // histEem_p->GetYaxis()->SetRangeUser(0,40000);
  // histEem_e->Fit("gaus");
  histEem_e->Draw();
  histEem_mu->DrawNormalized("same",histEdepo_e->Integral());
  histEem_had->DrawNormalized("same",histEdepo_e->Integral());
  // histEem_pi->Draw("same");
  // histEem_k->Draw("same");
  // histEem_p->Draw("same");
  // gPad->BuildLegend(0.7,0.8,0.95,0.95);
  gPad->BuildLegend(0.65,0.75,0.95,0.95);
  c1->cd(2);
  histEhad_e->SetLineColor(kRed);
  histEhad_mu->SetLineColor(kViolet);
  // histEem_pi->SetLineColor(kGreen);
  // histEem_k->SetLineColor(kViolet );

  // histEhad_mu->Scale(1.0/histEhad_mu->Integral());
  // histEhad_e->Scale(1.0/histEhad_e->Integral());
  // histEhad_had->Scale(1.0/histEhad_had->Integral());
  
  // histEhad_e->GetYaxis()->SetRangeUser(0,1);
  // histEhad_mu->GetYaxis()->SetRangeUser(0,1);
  // histEhad_had->GetYaxis()->SetRangeUser(0,1);
  // histEhad_pi->GetYaxis()->SetRangeUser(0,35000);
  // histEhad_k->GetYaxis()->SetRangeUser(0,35000);
  // histEhad_p->GetYaxis()->SetRangeUser(0,35000);

  histEhad_e->Draw();
  histEhad_mu->DrawNormalized("same",histEdepo_e->Integral());
  histEhad_had->DrawNormalized("same",histEdepo_e->Integral());

  // histEhad_pi->Draw("same");
  // histEhad_k->Draw("same");
  // histEhad_p->Draw("same");
  gPad->BuildLegend(0.65,0.75,0.95,0.95);
  c1->cd(3);
  histEdepo_e->SetLineColor(kRed);
  histEdepo_mu->SetLineColor(kViolet);
  // histEem_pi->SetLineColor(kGreen);
  // histEem_k->SetLineColor(kViolet );

  // histEdepo_mu->Scale(1.0/histEdepo_mu->Integral());
  // histEdepo_e->Scale(1.0/histEdepo_e->Integral());
  // histEdepo_had->Scale(1.0/histEdepo_had->Integral());

  // histEdepo_e->GetYaxis()->SetRangeUser(0,1);
  // histEdepo_mu->GetYaxis()->SetRangeUser(0,1);
  // histEdepo_had->GetYaxis()->SetRangeUser(0,1);

  // histEdepo_pi->GetYaxis()->SetRangeUser(0,25000);
  // histEdepo_k->GetYaxis()->SetRangeUser(0,25000);
  // histEdepo_p->GetYaxis()->SetRangeUser(0,25000);
  
  histEdepo_e->Draw();
  histEdepo_mu->DrawNormalized("same",histEdepo_e->Integral());
  histEdepo_had->DrawNormalized("same",histEdepo_e->Integral());
  // histEdepo_pi->Draw("same");
  // histEdepo_k->Draw("same");
  // histEdepo_p->Draw("same");
  gPad->BuildLegend(0.6,0.75,0.98,0.95);
  // c1->Divide(3,5);
  // c1->cd(1);
  // histEem_mu->Draw();
  // c1->cd(2);
  // histEhad_mu->Draw();
  // c1->cd(3);
  // histEdepo_mu->Draw();
  // c1->cd(4);
  // histEem_e->Draw();
  // c1->cd(5);
  // histEhad_e->Draw();
  // c1->cd(6);
  // histEdepo_e->Draw();
  // c1->cd(7);
  // histEem_pi->Draw();
  // c1->cd(8);
  // histEhad_pi->Draw();
  // c1->cd(9);
  // histEdepo_pi->Draw();
  // c1->cd(10);
  // histEem_k->Draw();
  // c1->cd(11);
  // histEhad_k->Draw();
  // c1->cd(12);
  // histEdepo_k->Draw();
  // c1->cd(13);
  // histEem_p->Draw();
  // c1->cd(14);
  // histEhad_p->Draw();
  // c1->cd(15);
  // histEdepo_p->Draw();
  c1->Print("fig/B_e1e1pi/Eem_Ehad.pdf"); 
  c1->Print("fig/B_e1e1pi/Eem_Ehad.png"); 
  delete c1;

  TCanvas *c2=new TCanvas("c2","c2",700,400);
  TGraph *Eem_reso = new TGraph(20,xE_em,xEem);
  TGraph *Ehad_reso = new TGraph(20,xE_had,xEhad);
  Eem_reso->GetYaxis()->SetTitle("#sigma_E/E for ecal");
  Eem_reso->GetXaxis()->SetTitle("E");
  Eem_reso->GetYaxis()->SetRangeUser(0,0.7);

  Ehad_reso->GetYaxis()->SetTitle("#sigma_E/E for hcal");
  Ehad_reso->GetXaxis()->SetTitle("E");
  Ehad_reso->GetYaxis()->SetRangeUser(0,0.7);

  c2->Divide(2,1);
  c2->cd(1);
  Eem_reso->Draw("APL");
  c2->cd(2);
  Ehad_reso->Draw("APL");
  c2->Print("fig/B_e1e1pi/Edepo_reso.png");
  c2->Print("fig/B_e1e1pi/Edepo_reso.pdf");
}


void Signal_Mass(ExRootTreeReader *treeReader){
  
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchPFC = treeReader->UseBranch("ParticleFlowCandidate");
    TH1F *hist_sig_mass = new TH1F("sig inv_mass", "sig inv_mass", 100,5,5.6);

    
    Long64_t numberOfEntries = treeReader->GetEntries();
    for(Long64_t entry = 0; entry <numberOfEntries ; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      if(branchJet->GetEntriesFast()!= 2) continue;
      Double_t Mass[2]={0.0,0.0};
      for(Int_t i=0;i<2;i++){
        Jet* jet=(Jet*) branchJet->At(i);
        cout<<"  n = "<<jet->Constituents.GetEntriesFast()<<endl;
      //  Mass[i]=Inv_Mass(jet);
      }
      if(Mass[0]>Mass[1]){
        hist_sig_mass->Fill(Mass[0]);
      }
      else{
        hist_sig_mass->Fill(Mass[1]);
      }
     
    }  
    TCanvas *c1=new TCanvas("c1","c1",700,700);
    hist_sig_mass->Draw();
    c1->Print("fig/B_e1e1pi/hist_mass.png"); 
    c1->Print("fig/B_e1e1pi/hist_mass.pdf"); 
    c1->Clear();
    delete c1;
}

void Bg_ana(ExRootTreeReader *treeReadersig,ExRootTreeReader *treeReaderbg){
  
    TClonesArray *branchJet_bg  = treeReaderbg ->UseBranch("Jet");
    TClonesArray *branchJet_sig = treeReadersig->UseBranch("Jet");
    TClonesArray *branchPFC_bg  = treeReaderbg ->UseBranch("ParticleFlowCandidate");
    TClonesArray *branchPFC_sig = treeReadersig->UseBranch("ParticleFlowCandidate");
    TH1F *hist_bg    = new TH1F("bg"   , "bg"   , 100,5,5.6);
    TH1F *hist_l_bg  = new TH1F("l_bg" , "l_bg" , 100,0,20 );
    TH1F *hist_sig   = new TH1F("sig"  , "sig"  , 100,5,5.6);
    TH1F *hist_l_sig = new TH1F("l_sig", "l_sig", 100,0,20 );

    double mass_b=5.27934;
    Long64_t numberOfEntries = treeReadersig->GetEntries();
    for(Long64_t entry = 0; entry <numberOfEntries ; ++entry)
    {
      // Load selected branches with data from specified event
      treeReadersig->ReadEntry(entry);
      if(branchJet_sig->GetEntriesFast()!= 2) continue;
      Double_t Mass[2]={0.0,0.0};
      Double_t mass=0.0;
      Double_t l_1=0.0,l_2=0.0,l=0.0;
      TLorentzVector P_pi,P_em,P_ep;
      TLorentzVector P_pi_1,P_em_1,P_ep_1;
      TLorentzVector P_pi_2,P_em_2,P_ep_2;
      for(Int_t i=0;i<2;i++){
        Jet* jet=(Jet*) branchJet_sig->At(i);
        if(i==0)  Mass[i]=Inv_Mass(jet,P_pi_1,P_em_1,P_ep_1,l_1);
        else      Mass[i]=Inv_Mass(jet,P_pi_2,P_em_2,P_ep_2,l_2);
      }
      if(abs(Mass[0]-mass_b)<abs(Mass[1]-mass_b)){
	mass=Mass[0];
        P_pi=P_pi_1;
        P_em=P_em_1;
        P_ep=P_ep_1;
        l=l_1;
      }
      else{
	mass=Mass[1];
        P_pi=P_pi_2;
        P_em=P_em_2;
        P_ep=P_ep_2;
        l=l_2;
      }
     if(l!=0) hist_l_sig->Fill(l);
     if(P_pi.Pt()<3 || P_em.Pt()<3 || P_ep.Pt()<3) continue;
     if(mass<5 || mass>5.6) continue;
     hist_sig->Fill(mass);
     
    }  

    for(Long64_t entry = 0; entry <treeReaderbg->GetEntries(); ++entry)
    {
      // Load selected branches with data from specified event
      treeReaderbg->ReadEntry(entry);
      if(branchJet_bg->GetEntriesFast()!= 2) continue;
      Double_t Mass[2]={0.0,0.0};
      Double_t mass=0.0;
      Double_t l_1=0.0,l_2=0.0,l=0.0;
      TLorentzVector P_pi,P_em,P_ep;
      TLorentzVector P_pi_1,P_em_1,P_ep_1;
      TLorentzVector P_pi_2,P_em_2,P_ep_2;
      for(Int_t i=0;i<2;i++){
        Jet* jet=(Jet*) branchJet_bg->At(i);
        if(i==0)  Mass[i]=Inv_Mass(jet,P_pi_1,P_em_1,P_ep_1,l_1);
        else      Mass[i]=Inv_Mass(jet,P_pi_2,P_em_2,P_ep_2,l_2);
      }
      if(abs(Mass[0]-mass_b)<abs(Mass[1]-mass_b)){
	mass=Mass[0];
        P_pi=P_pi_1;
        P_em=P_em_1;
        P_ep=P_ep_1;
        l=l_1;
      }
      else{
	mass=Mass[1];
        P_pi=P_pi_2;
        P_em=P_em_2;
        P_ep=P_ep_2;
        l=l_2;
      }
     if(l!=0) hist_l_bg->Fill(l);
     if(P_pi.Pt()<3 || P_em.Pt()<3 || P_ep.Pt()<3) continue;
     if(mass<5 || mass>5.6) continue;
     hist_bg->Fill(mass);

    }  
    TCanvas *c1=new TCanvas("c1","c1",700,700);
    gStyle->SetOptStat("im");
    hist_bg->Draw();
    c1->Print("fig/B_e1e1pi/hist_bg.png"); 
    c1->Print("fig/B_e1e1pi/hist_bg.pdf"); 
    hist_l_bg ->Draw();
    hist_l_sig->SetLineColor(kRed);
    hist_l_sig->DrawNormalized("same",hist_l_bg->Integral());
    gPad->BuildLegend(0.65,0.75,0.95,0.95);
    c1->Print("fig/B_e1e1pi/hist_l.png"); 
    c1->Print("fig/B_e1e1pi/hist_l.pdf"); 
    c1->Clear();
    delete c1;
}


double Inv_Mass(Jet *jet,TLorentzVector &P_pi,TLorentzVector &P_em,TLorentzVector &P_ep,double &l){
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
              if((canpi->X+canem->X)==0 || (canpi->Y+canem->Y)==0 || (canpi->Z+canem->Z)==0 || (canpi->T+canem->T)==0) continue;
              if((canep->X+canem->X)==0 || (canep->Y+canem->Y)==0 || (canep->Z+canem->Z)==0 || (canep->T+canem->T)==0) continue;
              if((canpi->X+canep->X)==0 || (canpi->Y+canep->Y)==0 || (canpi->Z+canep->Z)==0 || (canpi->T+canep->T)==0) continue;
              double l_pi_em = sqrt(pow((canpi->X-canem->X)/(canpi->X+canem->X),2)+pow((canpi->Y-canem->Y)/(canpi->Y+canem->Y),2)+pow((canpi->Z-canem->Z)/(canpi->Z+canem->Z),2)+pow((canpi->T-canem->T)/(canpi->T+canem->T),2));
              double l_ep_em = sqrt(pow((canep->X-canem->X)/(canep->X+canem->X),2)+pow((canep->Y-canem->Y)/(canep->Y+canem->Y),2)+pow((canep->Z-canem->Z)/(canep->Z+canem->Z),2)+pow((canep->T-canem->T)/(canep->T+canem->T),2));
              double l_pi_ep = sqrt(pow((canpi->X-canep->X)/(canpi->X+canep->X),2)+pow((canpi->Y-canep->Y)/(canpi->Y+canep->Y),2)+pow((canpi->Z-canep->Z)/(canpi->Z+canep->Z),2)+pow((canpi->T-canep->T)/(canpi->T+canep->T),2));
             l=sqrt(pow(l_pi_em,2)+pow(l_ep_em,2)+pow(l_pi_ep,2));
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


void Truth_Mass(ExRootTreeReader *treeReader){

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");

  TH1F *hist_truth_mass = new TH1F("B_truth_mass", "B_truth_mass", 100,5,5.6);

  TH1F *hist_l_sig = new TH1F("l_sig", "l_sig", 100,0,20 );

  GenParticle *particle;
  GenParticle *par;
  GenParticle *parpi;
  GenParticle *parep;
  GenParticle *parem;
  ParticleFlowCandidate *particleflowcandidate;
  ParticleFlowCandidate *canpi;
  ParticleFlowCandidate *canep;
  ParticleFlowCandidate *canem;
  Int_t B_nub = 0;
  Int_t num_ep = 0;
  Int_t num_em = 0;
  Int_t num_pi = 0;

  Long64_t numberOfEntries = treeReader->GetEntries();
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    treeReader->ReadEntry(entry);

    std::vector <Int_t> id_B;
    std::map <Int_t,Int_t> id_pi;
    std::map <Int_t,Int_t> id_em;
    std::map <Int_t,Int_t> id_ep;
    vector<Int_t>::iterator it_B;
    map<Int_t,Int_t>::iterator it_pi;
    map<Int_t,Int_t>::iterator it_em;
    map<Int_t,Int_t>::iterator it_ep;

    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      par=(GenParticle* )branchParticle->At(i);
      if(abs(par->PID)==521){
        id_B.push_back(i);
        B_nub++;
      }
    }
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);
      // particle = (GenParticle*) branchParticle->At(i);
      for(it_B = id_B.begin();it_B!=id_B.end();it_B++){
        if(particle->M1==*it_B && particleflowcandidate->truth_PID != 130){
          switch(particle->PID){
            case 211  : 
              id_pi.insert(pair<Int_t,Int_t>(distance(id_B.begin(),it_B),i));
              break;
            case -211  : 
              id_pi.insert(pair<Int_t,Int_t>(distance(id_B.begin(),it_B),i));
              break;
            case 11 :
              id_ep.insert(pair<Int_t,Int_t>(distance(id_B.begin(),it_B),i));
              break;
            case -11 :
              id_em.insert(pair<Int_t,Int_t>(distance(id_B.begin(),it_B),i));
              break;
          }
        }
      }
    }
    num_pi = num_pi+distance(id_pi.begin(),id_pi.end());
    num_ep = num_ep+distance(id_ep.begin(),id_ep.end());
    num_em = num_em+distance(id_em.begin(),id_em.end());
    // for(it_B = id_B.begin();it_B!=id_B.end();it_B++){
    //   for(it_pi=id_pi.begin();it_pi!=id_pi.end();it_pi++){
    //     if(it_pi->first ==distance(id_B.begin(),it_B)){}
    //   }
    // }
    double l=0.0;
    if(id_pi.empty() || id_ep.empty() || id_em.empty()) continue;
    for(it_pi=id_pi.begin();it_pi!=id_pi.end();it_pi++){
      canpi=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(it_pi->second);
      parpi=(GenParticle*) canpi->Particles.At(0);
      // parpi=(GenParticle*) branchParticle->At(it_pi->second);
      for(it_ep=id_ep.begin();it_ep!=id_ep.end();it_ep++){
        canep=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(it_ep->second);
        parep=(GenParticle*) canep->Particles.At(0);
        // parep=(GenParticle*) branchParticle->At(it_ep->second);
        for(it_em=id_em.begin();it_em!=id_em.end();it_em++){
          canem=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(it_em->second);
          parem=(GenParticle*) canem->Particles.At(0);
          // parem=(GenParticle*) branchParticle->At(it_em->second);
          if(it_pi->first == it_ep->first && it_ep->first == it_em->first){
            Double_t mass = (canpi->P4()+canep->P4()+canem->P4()).M();
            // cout<<"mass  ="<<mass<<endl;
            if(mass>5 && mass<5.6){
              hist_truth_mass->Fill(mass);
              double l_pi_em = sqrt(pow((parpi->X-parem->X)/(parpi->X+parem->X),2)+pow((parpi->Y-parem->Y)/(parpi->Y+parem->Y),2)+pow((parpi->Z-parem->Z)/(parpi->Z+parem->Z),2)+pow((parpi->T-parem->T)/(parpi->T+parem->T),2));
              double l_ep_em = sqrt(pow((parep->X-parem->X)/(parep->X+parem->X),2)+pow((parep->Y-parem->Y)/(parep->Y+parem->Y),2)+pow((parep->Z-parem->Z)/(parep->Z+parem->Z),2)+pow((parep->T-parem->T)/(parep->T+parem->T),2));
              double l_pi_ep = sqrt(pow((parpi->X-parep->X)/(parpi->X+parep->X),2)+pow((parpi->Y-parep->Y)/(parpi->Y+parep->Y),2)+pow((parpi->Z-parep->Z)/(parpi->Z+parep->Z),2)+pow((parpi->T-parep->T)/(parpi->T+parep->T),2));
             l=sqrt(pow(l_pi_em,2)+pow(l_ep_em,2)+pow(l_pi_ep,2));
             cout<<parpi->X<<endl;
             hist_l_sig->Fill(l);
            }
          }
        }   
      }
    }
    
    id_B.clear();
    id_pi.clear();
    id_ep.clear();
    id_em.clear();
  }
  cout<<"num_pi = "<<num_pi<<endl;
  cout<<"num_ep = "<<num_ep<<endl;
  cout<<"num_em = "<<num_em<<endl;
  cout<<"number of B is "<<B_nub<<endl;
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  hist_truth_mass->Draw();
  cout<<"integral = "<<hist_truth_mass->Integral()<<endl;
  c1->Print("fig/B_e1e1pi/hist_truth_mass.png"); 
  c1->Print("fig/B_e1e1pi/hist_truth_mass.pdf"); 
  hist_l_sig->Draw();
  c1->Print("fig/B_e1e1pi/hist_truth_l.png"); 
  c1->Clear();
  delete c1;
}
void TrackEfficiency(ExRootTreeReader *treeReader){
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");
  // TClonesArray *branchTrack = treeReader->UseBranch("Track");
  Long64_t numberOfEntries = treeReader->GetEntries();

  GenParticle *particle;
  // Track *track;
  ParticleFlowCandidate *particleflowcandidate;
  TH1F *histparticle_cos = new TH1F("truth cos", "truth cos ", 100,-1,1);
  TH1F *histeflowtrack_cos = new TH1F("eflowtrack cos", "eflowtrack cos", 100,-1,1);
  TH1F *histtrack_cos = new TH1F("track cos", "track cos", 100,-1,1);
  for(Int_t entry = 0; entry < 100000; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){
      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);
      /* if(particleflowcandidate->truth_PID == 22)   */histeflowtrack_cos->Fill(particleflowcandidate->CosTheta);
      // if(particle->PT < 0.5) continue;
    }
    // for(Long64_t i=0;i<branchTrack->GetEntriesFast();i++){
    //   track= (Track*) branchTrack->At(i);
    //   particle = (GenParticle*) track->Particle.GetObject();
    //   if(particle->PT < 0.5) continue;
    //   histtrack_cos->Fill(particle->P4().CosTheta());
    // }

    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      particle=(GenParticle* )branchParticle->At(i);
      if(particle->Charge!=0 && particle->Status == 1 ){
          // if(particle->PT < 0.5) continue;
          histparticle_cos->Fill(particle->P4().CosTheta());
      }
    }
  }
  TCanvas *c1=new TCanvas("c1","c1",600,600);
  // c1->Divide(2,1);
  // c1->cd(1);
  TH1D *hfrac = new TH1D("hfrac","h1/h2",100,-1,1);
  hfrac->Sumw2();
  hfrac->Divide(histeflowtrack_cos,histparticle_cos,1,1);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  histparticle_cos->SetLineColor(kRed);
  histparticle_cos->SetFillColor(kRed);
  histparticle_cos->GetYaxis()->SetRangeUser(0,4000);  
  // histparticle_cos->Draw("");
  histtrack_cos->SetLineColor(kGreen);
  histtrack_cos->SetFillColor(kGreen);
  histtrack_cos->GetYaxis()->SetRangeUser(0,4000);
  // histtrack_cos->Draw("same");
  // histeflowtrack_cos->GetYaxis()->SetRangeUser(0,4000);
  // histeflowtrack_cos->SetFillColor(kBlue);
  histeflowtrack_cos->Draw("same");
  // gPad->BuildLegend(0.7,0.8,0.95,0.95);


  // c1->Print("fig/B_e1e1pi/TrackEfficiency.png"); 
  // c1->Print("fig/B_e1e1pi/TrackEfficiency.pdf"); 
  // c1->cd(2);
  // hfrac->Draw();
  c1->Print("fig/B_e1e1pi/TrackEfficiency.png"); 
  
  delete c1;
}

void Performance_PID(ExRootTreeReader *treeReader){
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");
  Long64_t numberOfEntries = treeReader->GetEntries();

  GenParticle *particle;
  Track *track;
  ParticleFlowCandidate *particleflowcandidate;
  TH1D *hProb_Pi_pi = new TH1D("Prob_Pi for #pi", "Prob_Pi for #pi", 100,0,1.01);
  TH1D *hProb_Pi_k = new TH1D("Prob_Pi for k", "Prob_Pi for k", 100,0,1.01);
  TH1D *hProb_K_pi = new TH1D("Prob_K for #pi", "Prob_K for #pi", 100,0,1.01);
  TH1D *hProb_K_k = new TH1D("Prob_K for k", "Prob_K for k", 100,0,1.01);
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){
      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);
      if(abs(particleflowcandidate->truth_PID) == 211){
        hProb_K_pi->Fill(particleflowcandidate->Prob_K);
        hProb_Pi_pi->Fill(particleflowcandidate->Prob_Pi);
      }
      if(abs(particleflowcandidate->truth_PID) == 321){
        hProb_K_k->Fill(particleflowcandidate->Prob_K);
        hProb_Pi_k->Fill(particleflowcandidate->Prob_Pi);
      }

    }
  }
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  hProb_Pi_pi->SetLineColor(kRed);
  hProb_Pi_pi->SetLineWidth(3);
  hProb_Pi_pi->GetXaxis()->SetTitle("Prob_Pi");  
  hProb_Pi_pi->GetYaxis()->SetRangeUser(0,600000);  
  hProb_Pi_pi->Draw("");
  hProb_Pi_k->SetLineWidth(3);
  hProb_Pi_k->DrawNormalized("same",hProb_Pi_pi->Integral());
  gPad->BuildLegend(0.65,0.83,0.95,0.95);
  // c1->SetLogx();
  c1->Print("fig/B_e1e1pi/Prob_Pi.png"); 
  c1->Print("fig/B_e1e1pi/Prob_Pi.pdf"); 

  hProb_K_pi->SetLineColor(kRed);
  hProb_K_pi->SetLineWidth(3);
  hProb_K_pi->GetXaxis()->SetTitle("Prob_K");
  hProb_K_pi->GetYaxis()->SetRangeUser(0,600000);  
  hProb_K_pi->Draw("");
  hProb_K_k->SetLineWidth(3);
  hProb_K_k->DrawNormalized("same",hProb_K_pi->Integral()); 
  gPad->BuildLegend(0.65,0.83,0.95,0.95);
  c1->Print("fig/B_e1e1pi/Prob_K.png"); 
  c1->Print("fig/B_e1e1pi/Prob_K.pdf"); 
  delete c1;
  cout<<"Prob_pi for pi = "<<hProb_Pi_pi->Integral()<<endl;
  cout<<"Prob_k for pi = "<<hProb_K_pi->GetEntries()<<endl;
}

void Momentum(ExRootTreeReader *treeReader, ExRootTreeReader *treeReaderbg ){
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");
  TClonesArray *branchParticleFlowCandidate_bg = treeReaderbg->UseBranch("ParticleFlowCandidate");

  Long64_t numberOfEntries = treeReader->GetEntries();

  
  GenParticle *particle;
  TH1F *hmom_fB = new TH1F("sig", "sig", 100,0,10);
  TH1F *hmom_nfB = new TH1F("bg", "bg", 100,0,10);

  for(Long64_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::vector <Int_t> id_B;
    vector<Int_t>::iterator it_B;
    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      particle=(GenParticle* )branchParticle->At(i);
      if(abs(particle->PID)==521){
        id_B.push_back(i);
      }
    }
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){
      ParticleFlowCandidate *particleflowcandidate =(ParticleFlowCandidate*)branchParticleFlowCandidate->At(i);
      particle= (GenParticle*) particleflowcandidate->Particles.At(0);
      for(it_B = id_B.begin();it_B!=id_B.end();it_B++){
        if(particle->M1 == *it_B){
          hmom_fB->Fill(particleflowcandidate->P4().Pt());
        }
      }
    }
  }
  for(Long64_t entry = 0; entry < treeReaderbg->GetEntries(); ++entry)
  {
    // Load selected branches with data from specified event
    treeReaderbg->ReadEntry(entry);
    for(Long64_t i=0;i<branchParticleFlowCandidate_bg->GetEntriesFast();i++){
      ParticleFlowCandidate *particleflowcandidate =(ParticleFlowCandidate*)branchParticleFlowCandidate_bg->At(i);
        if(abs(particle->PID) == 11 || abs(particle->PID) == 211 ){
          hmom_nfB->Fill(particleflowcandidate->P4().Pt());
        }
      }
    }

  TCanvas *c1=new TCanvas("c1","c1",600,600);
  hmom_fB->SetLineColor(kRed);
  hmom_nfB->Draw();
  hmom_fB->DrawNormalized("same",hmom_nfB->Integral());
  gPad->BuildLegend(0.65,0.75,0.98,0.95);
  c1->Print("fig/B_e1e1pi/hmom.png"); 
  c1->Print("fig/B_e1e1pi/hmom.pdf"); 


  delete c1;
}
void testcos(ExRootTreeReader *treeReader){
 TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");
  Long64_t numberOfEntries = treeReader->GetEntries();

  GenParticle *particle;
  ParticleFlowCandidate *pfc;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    treeReader->ReadEntry(entry);
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){
      ParticleFlowCandidate *pfc =(ParticleFlowCandidate*)branchParticleFlowCandidate->At(i);
      TLorentzVector P4;
      P4.SetPtEtaPhiM(pfc->PT,pfc->Eta,pfc->Phi,pfc->Mass);
      cout<<P4.CosTheta()-pfc->CosTheta<<endl;

    }


  }  





}
