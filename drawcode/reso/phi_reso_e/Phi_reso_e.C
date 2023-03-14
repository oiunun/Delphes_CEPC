#define phi_reso_e_cxx
#include "phi_reso_e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void phi_reso_e::Loop()
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

  // Book histograms
//   TH2F *sigmaphi = new TH2F("small","#sigma#phi", 100, 0,80,100,0,0.02);
  TH1F *t_gaus1[30];       //each bin histogram
  for(Int_t i=0;i<30;i++){ char t_hname[20];
    sprintf(t_hname,"hg1%02d",i);
    t_gaus1[i]=new TH1F(t_hname,"",200,-0.05,0.05);}
//   TH1F *t_gaus2[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char t_hname[20];
//     sprintf(t_hname,"hg2%02d",i);
//     t_gaus2[i]=new TH1F(t_hname,"",200,0,0.02);}
//   TH1F *t_gaus3[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char t_hname[20];
//     sprintf(t_hname,"hg3%02d",i);
//     t_gaus3[i]=new TH1F(t_hname,"",200,0,0.02);}
//   TH1F *t_gaus4[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char t_hname[20];
//     sprintf(t_hname,"hg4%02d",i);
//     t_gaus4[i]=new TH1F(t_hname,"",200,0,0.02);}
//   TH1F *t_gaus5[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char t_hname[20];
//     sprintf(t_hname,"hg5%02d",i);
//     t_gaus5[i]=new TH1F(t_hname,"",200,0,0.02);}

  double t_aver1[30];
  double t_error1[30];
//   double t_aver2[30];
//   double t_error2[30];
//   double t_aver3[30];
//   double t_error3[30];
//   double t_aver4[30];
//   double t_error4[30];
//   double t_aver5[30];
//   double t_error5[30];

  // Loop over all events
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      for ( int j =0;  j<Photon_size; j++){
     // cout<<jentry<<" "<<j<<" "<<Track_ErrorD0[j]<<endl;
      Float_t var = Photon_Truth_Phi[j]-Photon_Phi[j];
      //  Float_t cut_var = abs(Track_CosTheta[j]);
       if(Photon_E[j]>1.) {
      //   if(var!=0. && cut_var>=0. && cut_var<=0.2 ) { 
        for ( int jp=0;jp<30;jp++){ if(Photon_E[j]>3.*jp && Photon_E[j]<=3.*jp+3.)
          {t_gaus1[jp]->Fill(var); } }/* } */
        
      //   if(var!=0. && cut_var>0.2 && cut_var <=0.4) {
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {t_gaus2[jp]->Fill(var); } }}
        
      //   if(var!=0. && cut_var>0.4 && cut_var <=0.6) { 
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {t_gaus3[jp]->Fill(var); } }}
        
      //   if(var!=0. && cut_var>0.6 && cut_var <=0.854) { 
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {t_gaus4[jp]->Fill(var); } }}
        
      //   if(var!=0. && cut_var>0.854 && cut_var <=0.98) { 
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {t_gaus5[jp]->Fill(var); } }}
        
     } } 
  }

  
  for(Int_t i=0;i<30;i++){
    /*fit[i]->SetParameters(500,gaus[i]->GetMean(),gaus[i]->GetRMS());
    gaus[i]->Fit(fit[i],"r");
    Double_t par[3];
    fit[i]->GetParameters(par);*/
    t_aver1[i]=t_gaus1[i]->GetMean();
    t_error1[i]=t_gaus1[i]->GetStdDev();
   //  t_aver2[i]=t_gaus2[i]->GetMean();
   //  t_error2[i]=t_gaus2[i]->GetStdDev();
   //  t_aver3[i]=t_gaus3[i]->GetMean();
   //  t_error3[i]=t_gaus3[i]->GetStdDev();
   //  t_aver4[i]=t_gaus4[i]->GetMean();
   //  t_error4[i]=t_gaus4[i]->GetStdDev();
   //  t_aver5[i]=t_gaus5[i]->GetMean();
   //  t_error5[i]=t_gaus5[i]->GetStdDev();
    delete t_gaus1[i];
   //  delete t_gaus2[i];
   //  delete t_gaus3[i];
   //  delete t_gaus4[i];
   //  delete t_gaus5[i];
   }

  TH1F *gaus1[30];       //each bin histogram
  for(Int_t i=0;i<30;i++){ char hname[20];
    sprintf(hname,"hg1%02d",i);
    gaus1[i]=new TH1F(hname,"",200,t_aver1[i]-5*t_error1[i],t_aver1[i]+5*t_error1[i]);}
  TH1F *pth1[30];        //each phi bin histogram
  for(Int_t i=0;i<30;i++){ char pthname[20];
    sprintf(pthname,"hpt1%02d",i);
    pth1[i]=new TH1F(pthname,"",100,3.*i,3.*i+3.);}
//   TH1F *gaus2[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char hname[20];
//     sprintf(hname,"hg2%02d",i);
//     gaus2[i]=new TH1F(hname,"",200,t_aver2[i]-5*t_error2[i],t_aver2[i]+5*t_error2[i]);}
//   TH1F *pth2[30];        //each pt bin histogram
//   for(Int_t i=0;i<30;i++){ char pthname[20];
//     sprintf(pthname,"hpt2%02d",i);
//     pth2[i]=new TH1F(pthname,"",100,3.*i,3.*i+3.);}
//   TH1F *gaus3[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char hname[20];
//     sprintf(hname,"hg3%02d",i);
//     gaus3[i]=new TH1F(hname,"",200,t_aver3[i]-5*t_error3[i],t_aver3[i]+5*t_error3[i]);}
//   TH1F *pth3[30];        //each pt bin histogram
//   for(Int_t i=0;i<30;i++){ char pthname[20];
//     sprintf(pthname,"hpt3%02d",i);
//     pth3[i]=new TH1F(pthname,"",100,3.*i,3.*i+3.);}
//   TH1F *gaus4[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char hname[20];
//     sprintf(hname,"hg4%02d",i);
//     gaus4[i]=new TH1F(hname,"",200,t_aver4[i]-5*t_error4[i],t_aver4[i]+5*t_error4[i]);}
//   TH1F *pth4[30];        //each pt bin histogram
//   for(Int_t i=0;i<30;i++){ char pthname[20];
//     sprintf(pthname,"hpt4%02d",i);
//     pth4[i]=new TH1F(pthname,"",100,3.*i,3.*i+3.);}
//    TH1F *gaus5[30];       //each bin histogram
//   for(Int_t i=0;i<30;i++){ char hname[20];
//     sprintf(hname,"hg5%02d",i);
//     gaus5[i]=new TH1F(hname,"",200,t_aver5[i]-5*t_error5[i],t_aver5[i]+5*t_error5[i]);}
//   TH1F *pth5[30];        //each pt bin histogram
//   for(Int_t i=0;i<30;i++){ char pthname[20];
//     sprintf(pthname,"hpt5%02d",i);
//     pth5[i]=new TH1F(pthname,"",100,3.*i,3.*i+3.);}

  double aver1[30];
  double error1[30];
  double errpt1[30];
  double pt1[30];
//   double aver2[30];
//   double error2[30];
//   double errpt2[30];
//   double pt2[30];
//   double aver3[30];
//   double error3[30];
//   double errpt3[30];
//   double pt3[30];
//   double aver4[30];
//   double error4[30];
//   double errpt4[30];
//   double pt4[30];
//   double aver5[30];
//   double error5[30];
//   double errpt5[30];
//   double pt5[30];


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      for ( int j =0;  j<Photon_size; j++){
     // cout<<jentry<<" "<<j<<" "<<Track_ErrorD0[j]<<endl;
      Float_t var = Photon_Truth_Phi[j]-Photon_Phi[j];
      //  Float_t cut_var = abs(Track_CosTheta[j]);
       if(Photon_E[j]>1.) {
      //   /* if(var!=0. && cut_var>=0. && cut_var<=0.2 ) { */ sigmaPT->Fill(Photon_E[j],var);
        for ( int jp=0;jp<30;jp++){ if(Photon_E[j]>3.*jp && Photon_E[j]<=3.*jp+3.)
          {gaus1[jp]->Fill(var); pth1[jp]->Fill(Photon_E[j]);} }/* } */
        
      //   if(var!=0. && cut_var>0.2 && cut_var <=0.4) { sigmaPT->Fill(Track_PT[j],var);
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {gaus2[jp]->Fill(var); pth2[jp]->Fill(Track_PT[j]);} }}
        
      //   if(var!=0. && cut_var>0.4 && cut_var <=0.6) { sigmaPT->Fill(Track_PT[j],var);
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {gaus3[jp]->Fill(var); pth3[jp]->Fill(Track_PT[j]);} }}
        
      //   if(var!=0. && cut_var>0.6 && cut_var <=0.854) { sigmaPT->Fill(Track_PT[j],var);
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {gaus4[jp]->Fill(var); pth4[jp]->Fill(Track_PT[j]);} }}
        
      //   if(var!=0. && cut_var>0.854 && cut_var <=0.98) { sigmaPT->Fill(Track_PT[j],var);
      //   for ( int jp=0;jp<30;jp++){ if(Track_PT[j]>3.*jp && Track_PT[j]<=3.*jp+3.)
      //     {gaus5[jp]->Fill(var); pth5[jp]->Fill(Track_PT[j]);} }}
        
     }}
   }

    for(Int_t i=0;i<30;i++){
    /*fit[i]->SetParameters(500,gaus[i]->GetMean(),gaus[i]->GetRMS());
    gaus[i]->Fit(fit[i],"r");
    Double_t par[3];
    fit[i]->GetParameters(par);*/
    pt1[i]=pth1[i]->GetMean();
    error1[i]=gaus1[i]->GetStdDev();
   //   pt2[i]=pth2[i]->GetMean();
   //  errpt2[i]=0.;
   //  aver2[i]=gaus2[i]->GetMean();
   //  error2[i]=gaus2[i]->GetMeanError();
   //   pt3[i]=pth3[i]->GetMean();
   //  errpt3[i]=0.;
   //  aver3[i]=gaus3[i]->GetMean();
   //  error3[i]=gaus3[i]->GetMeanError();
   //   pt4[i]=pth4[i]->GetMean();
   //  errpt4[i]=0.;
   //  aver4[i]=gaus4[i]->GetMean();
   //  error4[i]=gaus4[i]->GetMeanError();
   //   pt5[i]=pth5[i]->GetMean();
   //  errpt5[i]=0.;
   //  aver5[i]=gaus5[i]->GetMean();
   //  error5[i]=gaus5[i]->GetMeanError();
   //  std::cout<<gaus4[i]->GetEntries()<<" "<<aver4[i]<<" "<<gaus5[i]->GetEntries()<<" "<<aver5[i]<<endl;
   }
  
  // Show resulting histograms
  TGraph *g[5];
  g[0]=new TGraph(30,pt1,error1);
//   g[1]=new TGraphErrors(30,pt2,aver2,errpt2,error2);
//   g[2]=new TGraphErrors(24,pt3,aver3,errpt3,error3);
//   g[3]=new TGraphErrors(16,pt4,aver4,errpt4,error4);
//   g[4]=new TGraphErrors(12,pt5,aver5,errpt5,error5);

  TLegend *le=new TLegend(0.55, 0.6, 0.85, 0.85);
  int mak[] = {20, 21, 22, 23, 24, 25, 30};
  int line[] = {1, 1, 1, 1, 4, 4, 4};  
  int col[] = {kAzure+2, kOrange+1, kGreen+2, kMagenta+2, kOrange+1, kGreen+2, kMagenta+2};  
  for(int i=0;i<1;i++){
   //  char tit[200];
   //  if(0==i) sprintf(tit, "0<|cos#theta|<0.2");
   //  else if(1==i) sprintf(tit, "0.2<|cos#theta|<0.4");
   //  else if(2==i) sprintf(tit, "0.4<|cos#theta|<0.6");
   //  else if(3==i) sprintf(tit, "0.6<|cos#theta|<0.854");
   //  else if(4==i) sprintf(tit, "0.854<|cos#theta|<0.98");
    
    g[i]->SetMarkerStyle(mak[i]);
    g[i]->SetMarkerColor(col[i]);
    g[i]->SetMarkerSize(0.7);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineStyle(line[i]);
    g[i]->SetLineWidth(1);
   //  le->AddEntry(g[i],tit,"pl");
    }
  g[0]->GetYaxis()->SetTitle("#sigma#phi(rad)");
  g[0]->GetXaxis()->SetTitle("E(GeV)");
  g[0]->GetYaxis()->SetRangeUser(2e-3,5e-3);
  TCanvas *c1=new TCanvas("c1","c1",1300,700);
//   c1->SetLogy();
  g[0]->Draw("ALP");
//   for(int i=1;i<5;i++) g[i]->Draw("LP same");
//   le->Draw("same");
  c1->Print("../fig/Phi_reso_e.pdf");
  c1->Print("../fig/Phi_reso_e.png");
  // TCanvas*c2=new TCanvas("c2","c2");
  // sigmaPT->Draw();
  TFile *f=new TFile("../fig/root/his_Phi_e.root","recreate");
  for(int i = 0 ;i<30;i++){
    gaus1[i]->Write();
   //  gaus2[i]->Write();
   //  gaus3[i]->Write();
   //  gaus4[i]->Write();
   //  gaus5[i]->Write();
  }
  f->Close();
  return 0;


}
void Phi_reso_e(){
  phi_reso_e t;
  t.Loop();
  return 0;
}