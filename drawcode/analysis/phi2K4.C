/*
root -l examples/draw_invar_mass.C'("delphes_output.root")'
*/
#ifdef __CLING__
R__LOAD_LIBRARY(../../libDelphes.so)
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
#include "TMath.h"
#endif

//------------------------------------------------------------------------------
void Eem_Ehad(ExRootTreeReader *treeReader);
void Ana(ExRootTreeReader *treeReader,ExRootTreeReader *treeReaderbg);
void Sig_Bg_Mass(ExRootTreeReader *treeReader,ExRootTreeReader *treeReaderbg);
void Vphi2(Jet* jet,TLorentzVector &P_phi1,TLorentzVector &P_phi2);
void Truth_Mass(ExRootTreeReader *treeReader);
void TrackEfficiency(ExRootTreeReader *treeReader);
void Performance_PID(ExRootTreeReader *treeReader);
void Momentum(ExRootTreeReader *treeReader,ExRootTreeReader *treeReaderbg);

void phi2K4()
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  TChain chainbg("Delphes");
  chain.Add("../../../rootfile/phi2K4.root");
  chainbg.Add("../../../rootfile/mal/output_pythia6_bg/bb_*.root");
 // chain.Add("../../../rootfile/mal/output_pythia6/bb_*.root");
  // chain.Add("../../rootfile/qqh_X_1_9.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  ExRootTreeReader *treeReaderbg = new ExRootTreeReader(&chainbg);

  // Book histograms

  // Eem_Ehad(treeReader);
  Ana(treeReader,treeReaderbg);
  // Sig_Bg_Mass(treeReader,treeReaderbg);
 //  Truth_Mass(treeReader);
  // TrackEfficiency(treeReaderbg);
  // Performance_PID(treeReader);
//   Momentum(treeReader,treeReaderbg);

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


void Ana(ExRootTreeReader *treeReader,ExRootTreeReader *treeReaderbg){
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchParticleFlowCandidatebg = treeReaderbg->UseBranch("ParticleFlowCandidate");
  TClonesArray *branchJetbg = treeReaderbg->UseBranch("Jet");
  Long64_t numberOfEntries = treeReader->GetEntries();

  double mass_bs=5.36692; 
  TH1F *hist_sig_mass = new TH1F("sig inv_mass", "sig inv_mass",20,5,5.7);
  TH1F *hist_bg_mass = new TH1F("bg  inv_mass", "bg  inv_mass" ,20,5,5.7);
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    if(branchJet->GetEntriesFast()!= 2) continue;
    Double_t Mass[2]={0.0,0.0};
    Double_t mass=0.0;
    TLorentzVector P_phi1,P_phi2;
    TLorentzVector P_phi1_1,P_phi2_1;
    TLorentzVector P_phi1_2,P_phi2_2;
    for(Int_t i=0;i<2;i++){
      Jet* jet=(Jet*) branchJet->At(i);
      if(i==0){
        Vphi2(jet,P_phi1_1,P_phi2_1);
        Mass[i]=(P_phi1_1+P_phi2_1).M();
      }
      else{
        Vphi2(jet,P_phi1_2,P_phi2_2);
        Mass[i]=(P_phi1_2+P_phi2_2).M();
      }
     
    }
    if(abs(Mass[0]-mass_bs)<abs(Mass[1]-mass_bs)){
      mass=Mass[0];
      P_phi1=P_phi1_1;
      P_phi2=P_phi2_1;
    }
    else{
      mass=Mass[1];
      P_phi1=P_phi1_2;
      P_phi2=P_phi2_2;
    }
    if(mass>5 && mass<5.7 && abs(P_phi1.M()-1.019461)<0.025 && abs(P_phi1.M()-1.019461)<0.025){
      hist_sig_mass->Fill(mass);
    }
   }

  for(Int_t entry = 0; entry < treeReaderbg->GetEntries(); ++entry)
  {
    treeReaderbg->ReadEntry(entry);
    if(branchJetbg->GetEntriesFast()!= 2) continue;
    Double_t Mass[2]={0.0,0.0};
    Double_t mass=0.0;
    TLorentzVector P_phi1,P_phi2;
    TLorentzVector P_phi1_1,P_phi2_1;
    TLorentzVector P_phi1_2,P_phi2_2;
    for(Int_t i=0;i<2;i++){
      Jet* jet=(Jet*) branchJetbg->At(i);
      if(i==0){
        Vphi2(jet,P_phi1_1,P_phi2_1);
        Mass[i]=(P_phi1_1+P_phi2_1).M();
      }
      else{
        Vphi2(jet,P_phi1_2,P_phi2_2);
        Mass[i]=(P_phi1_2+P_phi2_2).M();
      }
     
    }
    if(abs(Mass[0]-mass_bs)<abs(Mass[1]-mass_bs)){
      mass=Mass[0];
      P_phi1=P_phi1_1;
      P_phi2=P_phi2_1;
    }
    else{
      mass=Mass[1];
      P_phi1=P_phi1_2;
      P_phi2=P_phi2_2;
    }
    if(mass>5 && mass<5.7 && abs(P_phi1.M()-1.019461)<0.025 && abs(P_phi1.M()-1.019461)<0.025){
      hist_bg_mass->Fill(mass);
    }
   }


   cout<<treeReaderbg->GetEntries()<<endl;
   TCanvas *c1=new TCanvas("c1","c1",700,700);
   gStyle->SetOptStat("mi");
    hist_sig_mass->Draw();
    c1->Print("fig/phi2K4/hist_sig_mass.png");
    hist_sig_mass->SetFillColor(kRed);
    hist_sig_mass->SetLineColor(kRed);
    hist_bg_mass->SetFillColor(kBlue);
    hist_sig_mass->Scale(1.0/10000.0*(float)treeReaderbg->GetEntries()/1109976.0);
    THStack *hs = new THStack("hs","Sig+Bg inv_mass");
    hs->Add(hist_bg_mass);
    hs->Add(hist_sig_mass);
    hs->Draw("hist");
    gPad->BuildLegend(0.75,0.85,0.95,0.95);
    c1->Print("fig/phi2K4/hist_mass.png"); 
    c1->Print("fig/phi2K4/hist_mass.pdf"); 
    
    c1->Clear();
    delete c1;

}


void Sig_Bg_Mass(ExRootTreeReader *treeReader,ExRootTreeReader *treeReaderbg){
  
    TH1F *hist_sig_mass = new TH1F("sig inv_mass", "sig ", 100,5,5.6);
    TH1F *hist_bg_mass = new TH1F("bg inv_mass", "bg ", 100,5,5.6);
 //   Inv_Mass(treeReader,hist_sig_mass);
 //   Inv_Mass(treeReaderbg,hist_bg_mass);
    hist_sig_mass->SetFillColor(kRed);
    hist_sig_mass->SetLineColor(kRed);
    hist_bg_mass->SetFillColor(kBlue);
    THStack *hs = new THStack("hs","Sig+Bg inv_mass");
    hs->Add(hist_bg_mass);
    hs->Add(hist_sig_mass);
    TCanvas *c1=new TCanvas("c1","c1",700,700);
    hs->Draw();
    gPad->BuildLegend(0.7,0.75,0.95,0.9);
    c1->Print("fig/B_e1e1pi/hmass_sig_bg.png"); 
    c1->Print("fig/B_e1e1pi/hmass_sig_bg.pdf"); 
    c1->Clear();
    delete c1;
}


void Vphi2(Jet* jet,TLorentzVector &P_phi1,TLorentzVector &P_phi2){
  ParticleFlowCandidate *particleflowcandidate;
  
  int n_Jetp = jet->Constituents.GetEntriesFast();
  if(n_Jetp==0) return 0;
  std::vector < Int_t > id_Kp;
  std::vector < Int_t > id_Km;
  for(Long64_t i=0;i<n_Jetp;i++){
    particleflowcandidate = (ParticleFlowCandidate*) jet->Constituents.At(i);
    if(particleflowcandidate->PID==321/* && particleflowcandidate->P>3 */){
      id_Kp.push_back(i);
    }
    if(particleflowcandidate->PID==-321/* && particleflowcandidate->P>3 */){
      id_Km.push_back(i);
    }
  }

  Double_t inv_mass = 0;
  TLorentzVector P_Kp1,P_Kp2,P_Km1,P_Km2;
  int i_Kp=0,i_Km=0;    
  if(!id_Kp.empty() && !id_Km.empty()) {    
    vector<int>::iterator it_Kp;
    vector<int>::iterator it_Km;
    for(it_Kp=id_Kp.begin();it_Kp!=id_Kp.end();it_Kp++){
      ParticleFlowCandidate *PFC_Kp = (ParticleFlowCandidate*) jet->Constituents.At(*it_Kp);
      for(it_Km=id_Km.begin();it_Km!=id_Km.end();it_Km++){
        ParticleFlowCandidate *PFC_Km = (ParticleFlowCandidate*) jet->Constituents.At(*it_Km);
        Double_t mass = (PFC_Kp->P4()+PFC_Km->P4()).M();
        if(abs(mass-1.019461)<abs(inv_mass-1.019461)) {
          inv_mass = mass;
          i_Kp = *it_Kp;
          i_Km = *it_Km;
          P_Kp1 = PFC_Kp->P4();
          P_Km1 = PFC_Km->P4();
        }
      }
    }
   //id_Kp.erase(find(id_Kp.begin(),id_Kp.end(),i_Kp));
  // id_Km.erase(find(id_Km.begin(),id_Km.end(),i_Km));
    }
  inv_mass = 0;
  if(!id_Kp.empty() && !id_Km.empty()) {       
    vector<int>::iterator it_Kp;
    vector<int>::iterator it_Km;
    for(it_Kp=id_Kp.begin();it_Kp!=id_Kp.end();it_Kp++){
      if(*it_Kp==i_Kp) continue;
      ParticleFlowCandidate *PFC_Kp = (ParticleFlowCandidate*) jet->Constituents.At(*it_Kp);
      for(it_Km=id_Km.begin();it_Km!=id_Km.end();it_Km++){
        if(*it_Km==i_Km) continue;
        ParticleFlowCandidate *PFC_Km = (ParticleFlowCandidate*) jet->Constituents.At(*it_Km);
        Double_t mass = (PFC_Kp->P4()+PFC_Km->P4()).M();
          if(abs(mass-1.019461)<abs(inv_mass-1.019461)) {
            inv_mass = mass;
            P_Kp2 = PFC_Kp->P4();
            P_Km2 = PFC_Km->P4();
          }
        }
      }
    }
  if(P_Kp1.Pt()<1.5 || P_Km1.Pt()<1.5 || P_Kp2.Pt()<1.5 || P_Km2.Pt()<1.5 ) {
    id_Kp.clear();
    id_Km.clear();
    return 0;
  }
   P_phi1=P_Kp1+P_Km1;
   P_phi2=P_Kp2+P_Km2;

   id_Kp.clear();
   id_Km.clear();

   return 0;
}


void Truth_Mass(ExRootTreeReader *treeReader){

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate = treeReader->UseBranch("ParticleFlowCandidate");

  TH1F *hist_truth_mass = new TH1F("Bs_truth_mass", "Bs_truth_mass", 100,5,5.7);

  GenParticle *particle;
  GenParticle *par;
  ParticleFlowCandidate *particleflowcandidate;
  ParticleFlowCandidate *canem;
  Int_t B_nub = 0;
  Long64_t numberOfEntries = treeReader->GetEntries();
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    treeReader->ReadEntry(entry);

    Int_t id_B=0,id_phi1=0,id_phi2=0,id_Kp1=0,id_Km1=0,id_Kp2=0,id_Km2=0;;
    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      par=(GenParticle* )branchParticle->At(i);
      if(abs(par->PID)==531){
        id_B=i;
        B_nub++;
      }
    }
    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      par=(GenParticle* )branchParticle->At(i);
      if(par->M1==id_B){
        if(id_phi1==0) id_phi1=i;
        else id_phi2=i;
      }
    }
    if(id_B==0 || id_phi1==0 || id_phi2==0) continue;
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);
      if(particle->M1==id_phi1 && particleflowcandidate->truth_PID != 130){
        if(particle->PID == 321)  id_Kp1=i;        
        if(particle->PID == -321) id_Km1=i;        
       }
      if(particle->M1==id_phi2 && particleflowcandidate->truth_PID != 130){
        if(particle->PID == 321)  id_Kp2=i;        
        if(particle->PID == -321) id_Km2=i;        
       }
      }
    if(id_Kp1==0 || id_Km1==0 || id_Kp2==0 || id_Km2==0) continue;
    ParticleFlowCandidate *PFC_Kp1=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(id_Kp1);
    ParticleFlowCandidate *PFC_Km1=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(id_Km1);
    ParticleFlowCandidate *PFC_Kp2=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(id_Kp2);
    ParticleFlowCandidate *PFC_Km2=(ParticleFlowCandidate*) branchParticleFlowCandidate->At(id_Km2);
    GenParticle *Par_Kp1=(GenParticle*) PFC_Kp1->Particles.At(0);
    GenParticle *Par_Km1=(GenParticle*) PFC_Km1->Particles.At(0);
    GenParticle *Par_Kp2=(GenParticle*) PFC_Kp2->Particles.At(0);
    GenParticle *Par_Km2=(GenParticle*) PFC_Km2->Particles.At(0);
    double m = (PFC_Kp1->P4()+PFC_Km1->P4()+PFC_Kp2->P4()+PFC_Km2->P4()).M();
    if(m>5.7 || m<5) continue;
    hist_truth_mass->Fill(m);  
    }
  cout<<"number of B is "<<B_nub<<endl;
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  hist_truth_mass->Draw();
  cout<<"integral = "<<hist_truth_mass->Integral()<<endl;
  c1->Print("fig/phi2K4/hist_truthPID_mass.png"); 
  c1->Print("fig/phi2K4/hist_truthPID_mass.pdf"); 
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
  TClonesArray *branchParticle_bg = treeReaderbg->UseBranch("Particle");
  TClonesArray *branchParticleFlowCandidate_bg = treeReaderbg->UseBranch("ParticleFlowCandidate");

  Long64_t numberOfEntries = treeReader->GetEntries();

  ParticleFlowCandidate *particleflowcandidate;
  GenParticle *particle;
  GenParticle *par;
  TH1F *hmom_sig = new TH1F("sig", "P of K ", 100,0,10);
  TH1F *hmom_bg  = new TH1F("bg" , "P of K ", 100,0,10);

  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    Int_t id_B=0,id_phi1=0,id_phi2=0;;
    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      par=(GenParticle* )branchParticle->At(i);
      if(abs(par->PID)==531){
        id_B=i;
      }
    }
    for(Long64_t i=0;i<branchParticle->GetEntriesFast();i++){
      par=(GenParticle* )branchParticle->At(i);
      if(par->M1==id_B){
        if(id_phi1==0) id_phi1=i;
        else id_phi2=i;
      }
    }
    if(id_B==0 || id_phi1==0 || id_phi2==0) continue;
    for(Long64_t i=0;i<branchParticleFlowCandidate->GetEntriesFast();i++){

      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);
      if((particle->M1==id_phi1 || particle->M1==id_phi2) && particleflowcandidate->truth_PID != 130){
        hmom_sig->Fill(particleflowcandidate->P4().Pt());
       }
    }
    
  }
  for(Int_t entry = 0; entry < treeReaderbg->GetEntries(); ++entry){
    treeReaderbg->ReadEntry(entry);
    for(Long64_t i=0;i<branchParticleFlowCandidate_bg->GetEntriesFast();i++){
      particleflowcandidate= (ParticleFlowCandidate*) branchParticleFlowCandidate_bg->At(i);
      particle = (GenParticle*) particleflowcandidate->Particles.At(0);
      if(abs(particle->PID)==321) hmom_bg->Fill(particleflowcandidate->P4().Pt());
    }  
    
  }
  TCanvas *c1=new TCanvas("c1","c1",600,600);
  gStyle->SetOptTitle(kFALSE);
  hmom_sig->SetLineColor(kRed);
  hmom_bg->Draw();
  hmom_sig->DrawNormalized("same",hmom_bg->Integral());
  gPad->BuildLegend(0.65,0.75,0.98,0.95);
  c1->Print("fig/phi2K4/hmom.png"); 
  c1->Print("fig/phi2K4/hmom.pdf"); 


  delete c1;
}
