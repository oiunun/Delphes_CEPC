#define d0_reso_cos_cxx
#include "d0_reso_cos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void d0_reso_cos::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ss16.C
//      root> ss16 t
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
  TH2F *sigmaPT = new TH2F("small","#sigmaPt-Pt", 100, 0,80,100,0,0.01);
  TH1F *t_gaus1[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char t_hname[20];
    sprintf(t_hname,"hg1%02d",i);
    t_gaus1[i]=new TH1F(t_hname,"",200,0,0.01);}
  TH1F *t_gaus2[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char t_hname[20];
    sprintf(t_hname,"hg2%02d",i);
    t_gaus2[i]=new TH1F(t_hname,"",200,0,0.01);}
  TH1F *t_gaus3[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char t_hname[20];
    sprintf(t_hname,"hg3%02d",i);
    t_gaus3[i]=new TH1F(t_hname,"",200,0,0.01);}
  TH1F *t_gaus4[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char t_hname[20];
    sprintf(t_hname,"hg4%02d",i);
    t_gaus4[i]=new TH1F(t_hname,"",200,0,0.01);}
  TH1F *t_gaus5[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char t_hname[20];
    sprintf(t_hname,"hg5%02d",i);
    t_gaus5[i]=new TH1F(t_hname,"",200,0,0.01);}

  double t_aver1[26];
  double t_error1[26];
  double t_aver2[26];
  double t_error2[26];
  double t_aver3[26];
  double t_error3[26];
  double t_aver4[26];
  double t_error4[26];
  double t_aver5[26];
  double t_error5[26];

  // Loop over all events
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      for ( int j =0;  j<Track_size; j++){
     // cout<<jentry<<" "<<j<<" "<<Track_ErrorD0[j]<<endl;
       Float_t var=Track_ErrorD0[j]; 
       Float_t abs_cos=abs(Track_CosTheta[j]);
       Float_t cut_var=Track_PT[j];
      //  if(Track_PT[j]>1.) {
        if(var!=0. && cut_var>=1.0 && cut_var<=3.0 ) { 
        for ( int jp=0;jp<26;jp++){ if( abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {t_gaus1[jp]->Fill(var); } }}
        
        if(var!=0. && cut_var>3.0 && cut_var<=10.) {
        for ( int jp=0;jp<26;jp++){ if( abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {t_gaus2[jp]->Fill(var); } }}
        
        if(var!=0. && cut_var>20. && cut_var<=30.) { 
        for ( int jp=0;jp<26;jp++){ if( abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {t_gaus3[jp]->Fill(var); } }}
        
        if(var!=0. && cut_var>30. && cut_var<=50.) { 
        for ( int jp=0;jp<26;jp++){ if( abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {t_gaus4[jp]->Fill(var); } }}
        
        if(var!=0. && cut_var> 60.&& cut_var<=80.) { 
        for ( int jp=0;jp<26;jp++){ if( abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {t_gaus5[jp]->Fill(var); } }}
        
   //   }
   }
  }

  
  for(Int_t i=0;i<26;i++){
    /*fit[i]->SetParameters(500,gaus[i]->GetMean(),gaus[i]->GetRMS());
    gaus[i]->Fit(fit[i],"r");
    Double_t par[3];
    fit[i]->GetParameters(par);*/
    t_aver1[i]=t_gaus1[i]->GetMean();
    t_error1[i]=t_gaus1[i]->GetStdDev();
    t_aver2[i]=t_gaus2[i]->GetMean();
    t_error2[i]=t_gaus2[i]->GetStdDev();
    t_aver3[i]=t_gaus3[i]->GetMean();
    t_error3[i]=t_gaus3[i]->GetStdDev();
    t_aver4[i]=t_gaus4[i]->GetMean();
    t_error4[i]=t_gaus4[i]->GetStdDev();
    t_aver5[i]=t_gaus5[i]->GetMean();
    t_error5[i]=t_gaus5[i]->GetStdDev();
    delete t_gaus1[i];
    delete t_gaus2[i];
    delete t_gaus3[i];
    delete t_gaus4[i];
    delete t_gaus5[i];
   }

  TH1F *gaus1[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char hname[20];
    sprintf(hname,"hg1%02d",i);
    if(i<20) gaus1[i]=new TH1F(hname,"",200,t_aver1[i]-10*t_error1[i],t_aver1[i]+10*t_error1[i]);
    else if(19<i && i<25) gaus1[i]=new TH1F(hname,"",200,t_aver1[i]-21*t_error1[i],t_aver1[i]+21*t_error1[i]);
    else if(25==i) gaus1[i]=new TH1F(hname,"",200,t_aver1[i]-20*t_error1[i],t_aver1[i]+70*t_error1[i]);}
  TH1F *pth1[26];        //each pt bin histogram
  for(Int_t i=0;i<26;i++){ char pthname[20];
    sprintf(pthname,"hpt1%02d",i);
    pth1[i]=new TH1F(pthname,"",100,0.03846*i,0.03846*(i+1.0));}
  TH1F *gaus2[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char hname[20];
    sprintf(hname,"hg2%02d",i);
    gaus2[i]=new TH1F(hname,"",200,t_aver2[i]-10*t_error2[i],t_aver2[i]+10*t_error2[i]);}
  TH1F *pth2[26];        //each pt bin histogram
  for(Int_t i=0;i<26;i++){ char pthname[20];
    sprintf(pthname,"hpt2%02d",i);
    pth2[i]=new TH1F(pthname,"",100,0.03846*i,0.03846*(i+1.0));}
  TH1F *gaus3[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char hname[20];
    sprintf(hname,"hg3%02d",i);
    gaus3[i]=new TH1F(hname,"",200,t_aver3[i]-10*t_error3[i],t_aver3[i]+10*t_error3[i]);}
  TH1F *pth3[26];        //each pt bin histogram
  for(Int_t i=0;i<26;i++){ char pthname[20];
    sprintf(pthname,"hpt3%02d",i);
    pth3[i]=new TH1F(pthname,"",100,0.03846*i,0.03846*(i+1.0));}
  TH1F *gaus4[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char hname[20];
    sprintf(hname,"hg4%02d",i);
    gaus4[i]=new TH1F(hname,"",200,t_aver4[i]-10*t_error4[i],t_aver4[i]+10*t_error4[i]);}
  TH1F *pth4[26];        //each pt bin histogram
  for(Int_t i=0;i<26;i++){ char pthname[20];
    sprintf(pthname,"hpt4%02d",i);
    pth4[i]=new TH1F(pthname,"",100,0.03846*i,0.03846*(i+1.0));}
   TH1F *gaus5[26];       //each bin histogram
  for(Int_t i=0;i<26;i++){ char hname[20];
    sprintf(hname,"hg5%02d",i);
    gaus5[i]=new TH1F(hname,"",200,t_aver5[i]-10*t_error5[i],t_aver5[i]+10*t_error5[i]);}
  TH1F *pth5[26];        //each pt bin histogram
  for(Int_t i=0;i<26;i++){ char pthname[20];
    sprintf(pthname,"hpt5%02d",i);
    pth5[i]=new TH1F(pthname,"",100,0.03846*i,0.03846*(i+1.0));}

  double aver1[26];
  double error1[26];
  double errpt1[26];
  double pt1[26];
  double aver2[26];
  double error2[26];
  double errpt2[26];
  double pt2[26];
  double aver3[26];
  double error3[26];
  double errpt3[26];
  double pt3[26];
  double aver4[26];
  double error4[26];
  double errpt4[26];
  double pt4[26];
  double aver5[26];
  double error5[26];
  double errpt5[26];
  double pt5[26];


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      for ( int j =0;  j<Track_size; j++){
     // cout<<jentry<<" "<<j<<" "<<Track_ErrorD0[j]<<endl;
       Float_t var=Track_ErrorD0[j]; 
       Float_t abs_cos=abs(Track_CosTheta[j]);
       Float_t cut_var=Track_PT[j];

        if(var!=0. && cut_var> 1.0 && cut_var<=3.0 ) { sigmaPT->Fill(Track_PT[j],var);
        for ( int jp=0;jp<26;jp++){ if(abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {gaus1[jp]->Fill(var); pth1[jp]->Fill(abs_cos);} }}
        
        if(var!=0. && cut_var>=3.0 && cut_var<=10. ) { sigmaPT->Fill(Track_PT[j],var);
        for ( int jp=0;jp<26;jp++){ if(abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {gaus2[jp]->Fill(var); pth2[jp]->Fill(abs_cos);} }}
        
        if(var!=0. && cut_var>20. && cut_var<=30.) { sigmaPT->Fill(Track_PT[j],var);
        for ( int jp=0;jp<26;jp++){ if(abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {gaus3[jp]->Fill(var); pth3[jp]->Fill(abs_cos);} }}
        
        if(var!=0. && cut_var>30. && cut_var<=50.) { sigmaPT->Fill(Track_PT[j],var);
        for ( int jp=0;jp<26;jp++){ if(abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {gaus4[jp]->Fill(var); pth4[jp]->Fill(abs_cos);} }}
        
        if(var!=0. && cut_var> 60.&& cut_var<=80.) { sigmaPT->Fill(Track_PT[j],var);
        for ( int jp=0;jp<26;jp++){ if(abs_cos>0.03846*jp && abs_cos<=0.03846*(jp+1.0))
          {gaus5[jp]->Fill(var); pth5[jp]->Fill(abs_cos);} }}
        
     }
   }

    for(Int_t i=0;i<26;i++){
    /*fit[i]->SetParameters(500,gaus[i]->GetMean(),gaus[i]->GetRMS());
    gaus[i]->Fit(fit[i],"r");
    Double_t par[3];
    fit[i]->GetParameters(par);*/
    pt1[i]=pth1[i]->GetMean();
    errpt1[i]=0.;
    aver1[i]=gaus1[i]->GetMean();
    error1[i]=gaus1[i]->GetMeanError();
     pt2[i]=pth2[i]->GetMean();
    errpt2[i]=0.;
    aver2[i]=gaus2[i]->GetMean();
    error2[i]=gaus2[i]->GetMeanError();
     pt3[i]=pth3[i]->GetMean();
    errpt3[i]=0.;
    aver3[i]=gaus3[i]->GetMean();
    error3[i]=gaus3[i]->GetMeanError();
     pt4[i]=pth4[i]->GetMean();
    errpt4[i]=0.;
    aver4[i]=gaus4[i]->GetMean();
    error4[i]=gaus4[i]->GetMeanError();
     pt5[i]=pth5[i]->GetMean();
    errpt5[i]=0.;
    aver5[i]=gaus5[i]->GetMean();
    error5[i]=gaus5[i]->GetMeanError();
    std::cout<<gaus4[i]->GetEntries()<<" "<<aver4[i]<<" "<<error4[i]<<" "<<gaus5[i]->GetEntries()<<" "<<aver5[i]<<endl;
   }
  
  // Show resulting histograms
  TGraphErrors *g[5];
  g[0]=new TGraphErrors(26,pt1,aver1,errpt1,error1);
  g[1]=new TGraphErrors(26,pt2,aver2,errpt2,error2);
  g[2]=new TGraphErrors(26,pt3,aver3,errpt3,error3);
  g[3]=new TGraphErrors(24,pt4,aver4,errpt4,error4);
  g[4]=new TGraphErrors(17,pt5,aver5,errpt5,error5);

  TLegend *le=new TLegend(0.75, 0.15, 0.88, 0.3);
  int mak[] = {20, 21, 22, 23, 24, 25, 26};
  int line[] = {1, 1, 1, 1, 4, 4, 4};  
  int col[] = {kAzure+2, kOrange+1, kGreen+2, kMagenta+2, kOrange+1, kGreen+2, kMagenta+2};  
  for(int i=0;i<5;i++){
    char tit[200];
    if(0==i) sprintf(tit,"1.0<pt<3.0");
    else if(1==i) sprintf(tit, "3.0<pt<10.0");
    else if(2==i) sprintf(tit,  "20<pt<30");
    else if(3==i) sprintf(tit,  "30<pt<50");
    else if(4==i) sprintf(tit,  "50<pt<80");
    
    g[i]->SetMarkerStyle(mak[i]);
    g[i]->SetMarkerColor(col[i]);
    g[i]->SetMarkerSize(0.7);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineStyle(line[i]);
    g[i]->SetLineWidth(1);
    le->AddEntry(g[i],tit,"pl");
    }
  g[0]->GetYaxis()->SetTitle("#sigmaD0(cm)");
  g[0]->GetXaxis()->SetTitle("|cos#theta|");
  g[0]->GetYaxis()->SetRangeUser(0.001,0.02);
  TCanvas *c1=new TCanvas("c1","c1",1300,700);
  c1->SetLogy();
  g[0]->Draw("ALP");
  for(int i=1;i<5;i++) g[i]->Draw("LP same");
  le->Draw("same");
  c1->Print("../fig/d0_reso_cos.pdf");
  c1->Print("../fig/d0_reso_cos.png");
  // TCanvas*c2=new TCanvas("c2","c2");
  // sigmaPT->Draw();
  TFile *f=new TFile("../fig/root/his_d0_cos.root","recreate");
  for(int i = 0 ;i<26;i++){
    gaus1[i]->Write();
    gaus2[i]->Write();
    gaus3[i]->Write();
    gaus4[i]->Write();
    gaus5[i]->Write();
  }
  f->Close();



}
void D0_reso_cos(){
  d0_reso_cos t;
  t.Loop();
}
