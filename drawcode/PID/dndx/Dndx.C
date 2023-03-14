#define dndx_cxx
#include "dndx.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

Double_t Nclusters(Double_t begam, Int_t Opt);
void dndx::Loop()
{
//   In a ROOT session, you can do:
//      root> .L dndx.C
//      root> dndx t
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

   Long64_t nbytes = 0, nb = 0;
   TH1F* his_f = new TH1F("his","his",100,0,150);
   TFile *file=new TFile("e1e2p.root","recreate");
   TH1F *sigma1[25];
   for(Int_t i=0;i<25;i++){ 
      char t_hname[50];
      Float_t m0 = pow(10,-0.2+i*0.09);
      Float_t m1 = pow(10,-0.2+(i+1)*0.09);
      sprintf(t_hname,"%.3f<bg1<%.3f",m0,m1);
      sigma1[i]=new TH1F(t_hname,"#sigma' distribution",200,0,200);}
       TH1F *sigma2[25];
   for(Int_t i=0;i<25;i++){ 
      char t_hname[50];
      Float_t m0 = pow(10,-0.2+i*0.09);
      Float_t m1 = pow(10,-0.2+(i+1)*0.09);
      sprintf(t_hname,"%.3f<bg2<%.3f",m0,m1);
      sigma2[i]=new TH1F(t_hname,"#sigma' distribution",200,0,200);}
       TH1F *sigma3[25];
   for(Int_t i=0;i<25;i++){ 
      char t_hname[50];
      Float_t m0 = pow(10,-0.2+i*0.09);
      Float_t m1 = pow(10,-0.2+(i+1)*0.09);
      sprintf(t_hname,"%.3f<bg3<%.3f",m0,m1);
      sigma3[i]=new TH1F(t_hname,"#sigma' distribution",200,0,200);}
       TH1F *sigma4[25];
   for(Int_t i=0;i<25;i++){ 
      char t_hname[50];
      Float_t m0 = pow(10,-0.2+i*0.09);
      Float_t m1 = pow(10,-0.2+(i+1)*0.09);
      sprintf(t_hname,"%.3f<bg4<%.3f",m0,m1);
      sigma4[i]=new TH1F(t_hname,"#sigma' distribution",200,0,200);}
       TH1F *sigma5[25];
   for(Int_t i=0;i<25;i++){ 
      char t_hname[50];
      Float_t m0 = pow(10,-0.2+i*0.09);
      Float_t m1 = pow(10,-0.2+(i+1)*0.09);
      sprintf(t_hname,"%.3f<bg5<%.3f",m0,m1);
      sigma5[i]=new TH1F(t_hname,"#sigma' distribution",200,0,200);}


   TH1F *pth1[25];      
  for(Int_t i=0;i<25;i++){ char pthname[25];
    sprintf(pthname,"hpt1%02d",i);
    pth1[i]=new TH1F(pthname,"",100,pow(10,-0.2+i*0.09),pow(10,-0.2+(i+1)*0.09));}
       TH1F *pth2[25];      
  for(Int_t i=0;i<25;i++){ char pthname[25];
    sprintf(pthname,"hpt2%02d",i);
    pth2[i]=new TH1F(pthname,"",100,pow(10,-0.2+i*0.09),pow(10,-0.2+(i+1)*0.09));}
       TH1F *pth3[25];      
  for(Int_t i=0;i<25;i++){ char pthname[25];
    sprintf(pthname,"hpt3%02d",i);
    pth3[i]=new TH1F(pthname,"",100,pow(10,-0.2+i*0.09),pow(10,-0.2+(i+1)*0.09));}
       TH1F *pth4[25];      
  for(Int_t i=0;i<25;i++){ char pthname[25];
    sprintf(pthname,"hpt4%02d",i);
    pth4[i]=new TH1F(pthname,"",100,pow(10,-0.2+i*0.09),pow(10,-0.2+(i+1)*0.09));}
       TH1F *pth5[25];      
  for(Int_t i=0;i<25;i++){ char pthname[25];
    sprintf(pthname,"hpt5%02d",i);
    pth5[i]=new TH1F(pthname,"",100,pow(10,-0.2+i*0.09),pow(10,-0.2+(i+1)*0.09));}





   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for(Long64_t j=0;j<Track_size;j++){
         if(( Track_Nclusters[j]>0 && Track_Nclusters[j]<10000 && Track_L_DC[j]>0 && abs(Track_Truth_CosTheta[j])<0.98)&&(abs(Track_Truth_PID[j])==321 ||abs(Track_Truth_PID[j])==211 || abs(Track_Truth_PID[j])==2112)){
         Double_t x,y,eff,sigma,truth_ncl,bg;
         bg = Track_Truth_P[j]/Track_Mass[j];
         truth_ncl =Nclusters(bg,0)*Track_L_DC[j];   
         TH1F* h = new TH1F("h","h",100,truth_ncl-20*Track_Nclusters_err[j],truth_ncl+20*Track_Nclusters_err[j]);
         

         for (int i = 0; i < 1000;i++){
            x = gRandom->Gaus(truth_ncl,TMath::Sqrt(truth_ncl));
            eff = x*0.01*(-0.007309)/(Track_L_DC[j]*TMath::Sqrt(1-pow(Track_Truth_CosTheta[j],2)))+1.245497;
            y = gRandom->Gaus(eff,0.02); 
            h->Fill(x*y);
            
         }
         // h->Write();
         // cout<<"x*y = "<<truth_ncl <<endl;
         // cout<<"ncl "<<Track_Nclusters[j] <<endl;
         h->Fit("gaus","Q");//Q:Quiet mode (minimum printing)  0:拟合后不绘制直方图和拟合函数，但与选项"N"相反，它将拟合函数存储在函数直方图列表中
         TF1 *f =(TF1 *)h->GetFunction("gaus");
         sigma = f->GetParameter(2);

         Double_t eff_w,sigma_w;
         eff_w = truth_ncl*0.01*(-0.007309)/(Track_L_DC[j]*TMath::Sqrt(1-pow(Track_Truth_CosTheta[j],2)))+1.245497;
         sigma_w = eff_w*TMath::Sqrt(truth_ncl);

         Double_t cos =abs(Track_Truth_CosTheta[j]);
         for ( int jp=0;jp<25;jp++){ 
          if(bg>pow(10,-0.2+jp*0.09) && bg<pow(10,-0.2+(jp+1)*0.09)){
            if(0<cos<0.2) sigma1[jp]->Fill(sigma);pth1[jp]->Fill(bg);
            if(0.2<cos<0.4) sigma2[jp]->Fill(sigma);pth2[jp]->Fill(bg);
            if(0.4<cos<0.6) sigma3[jp]->Fill(sigma);pth3[jp]->Fill(bg);
            if(0.6<cos<0.8) sigma4[jp]->Fill(sigma);pth4[jp]->Fill(bg);
            if(0.8<cos<0.98) sigma5[jp]->Fill(sigma);pth5[jp]->Fill(bg);
            }
        }

         // his_f->Fill(sigma/sigma_w);
         delete f;
         delete h;
         }
      }
   }
   

  for(int i=0;i<25;i++){
  //  aver1[i]=sigma[i]->GetMean();
      sigma1[i]->Write();
      sigma2[i]->Write();
      sigma3[i]->Write();
      sigma4[i]->Write();
      sigma5[i]->Write();
  //  pt1[i]=pth1[i]->GetMean();
  //  cout<<"aver1 = "<<aver1[i]<<endl;
  //  cout<<"pth1 = "<<pt1[i]<<endl;
  }
  double aver1[25];
  double bg1[25];
  double aver2[25];
  double bg2[25];
  double aver3[25];
  double bg3[25];
  double aver4[25];
  double bg4[25];
  double aver5[25];
  double bg5[25];

  for(Int_t i=0;i<25;i++){
    /*fit[i]->SetParameters(500,gaus[i]->GetMean(),gaus[i]->GetRMS());
    gaus[i]->Fit(fit[i],"r");
    Double_t par[3];
    fit[i]->GetParameters(par);*/
    bg1[i]=pth1[i]->GetMean();
    aver1[i]=sigma1[i]->GetMean();
     bg2[i]=pth2[i]->GetMean();
    aver2[i]=sigma2[i]->GetMean();
     bg3[i]=pth3[i]->GetMean();
    aver3[i]=sigma3[i]->GetMean();
     bg4[i]=pth4[i]->GetMean();
    aver4[i]=sigma4[i]->GetMean();
     bg5[i]=pth5[i]->GetMean();
    aver5[i]=sigma5[i]->GetMean();
   }


  TGraph *g[5];
  g[0]=new TGraph(25,bg1,aver1);
  g[1]=new TGraph(25,bg2,aver2);
  g[2]=new TGraph(25,bg3,aver3);
  g[3]=new TGraph(25,bg4,aver4);
  g[4]=new TGraph(25,bg5,aver5);

  TLegend *le=new TLegend(0.75, 0.75, 0.88, 0.9);
  int mak[] = {20, 21, 22, 23, 24, 25, 26};
  int line[] = {1, 1, 1, 1, 4, 4, 4};  
  int col[] = {kAzure+2, kOrange+1, kGreen+2, kMagenta+2, kOrange+1, kGreen+2, kMagenta+2};  
  for(int i=0;i<5;i++){
    char tit[200];
    if(0==i) sprintf(tit,"0<cos#theta<0.2");
    else if(1==i) sprintf(tit, "0.2<cos#theta<0.4");
    else if(2==i) sprintf(tit,  "0.4<cos#theta<0.6");
    else if(3==i) sprintf(tit,  "0.6<cos#theta<0.8");
    else if(4==i) sprintf(tit,  "0.8<cos#theta<0.98");
    
    g[i]->SetMarkerStyle(mak[i]);
    g[i]->SetMarkerColor(col[i]);
    g[i]->SetMarkerSize(0.7);
    g[i]->SetLineColor(col[i]);
    g[i]->SetLineStyle(line[i]);
    g[i]->SetLineWidth(1);
    le->AddEntry(g[i],tit,"pl");
    }
  g[0]->GetYaxis()->SetTitle("#sigma_ncl");
  g[0]->GetXaxis()->SetTitle("#beta#gamma");
  g[0]->GetYaxis()->SetRangeUser(40,150);
  TCanvas *c1=new TCanvas("c1","c1",1300,700);
  // c1->SetLogy();
  g[1]->Draw("ALP");
  for(int i=2;i<5;i++) g[i]->Draw("LP same");
  le->Draw("same");
  c1->Print("sigma_ncl.pdf");
  c1->Print("sigma_ncl.png");

   file->Close();

  //  TGraph *g = new TGraph(25,pt1,aver1);
  //  TCanvas *c1=new TCanvas("c1","c1",700,700);
  //  g->SetMarkerStyle(20);

  //  g->GetYaxis()->SetTitle("#sigma'");
  //  g->GetXaxis()->SetTitle("bg");
  //  g->Draw();
  //  c1->SetLogx();
  //  c1->Print("sigma_dndx.png");
  //  c1->Print("sigma_dndx.pdf");
}
Double_t Nclusters(Double_t begam, Int_t Opt) {
	//
	// Opt = 0: He 90 - Isobutane 10
	//     = 1: pure He
	//     = 2: Argon 50 - Ethane 50
	//     = 3: pure Argon
	//
	//
	const Int_t Npt = 18;
	Double_t bg[Npt] = { 0.5, 0.8, 1., 2., 3., 4., 5., 8., 10.,
	12., 15., 20., 50., 100., 200., 500., 1000., 10000. };
	//
	// He 90 - Isobutane 10
	Double_t ncl_He_Iso[Npt] = { 42.94, 23.6,18.97,12.98,12.2,12.13,
	12.24,12.73,13.03,13.29,13.63,14.08,15.56,16.43,16.8,16.95,16.98, 16.98 };
	//
	// pure He
	Double_t ncl_He[Npt] = { 11.79,6.5,5.23,3.59,3.38,3.37,3.4,3.54,3.63,
				3.7,3.8,3.92,4.33,4.61,4.78,4.87,4.89, 4.89 };
	//
	// Argon 50 - Ethane 50
	Double_t ncl_Ar_Eth[Npt] = { 130.04,71.55,57.56,39.44,37.08,36.9,
	37.25,38.76,39.68,40.49,41.53,42.91,46.8,48.09,48.59,48.85,48.93,48.93 };
	//
	// pure Argon
	Double_t ncl_Ar[Npt] = { 88.69,48.93,39.41,27.09,25.51,25.43,25.69,
	26.78,27.44,28.02,28.77,29.78,32.67,33.75,34.24,34.57,34.68, 34.68 };
	//
	Double_t ncl[Npt];
	switch (Opt)
	{
	case 0: std::copy(ncl_He_Iso, ncl_He_Iso + Npt, ncl);	// He-Isobutane
		break;
	case 1: std::copy(ncl_He, ncl_He + Npt, ncl);		// pure He
		break;
	case 2: std::copy(ncl_Ar_Eth, ncl_Ar_Eth + Npt, ncl);	// Argon - Ethane
		break;
	case 3: std::copy(ncl_Ar, ncl_Ar + Npt, ncl);		// pure Argon
		break;
	}
	//
	Double_t interp = 0.0;
	TSpline3* sp3 = new TSpline3("sp3", bg, ncl, Npt);
	if (begam > bg[0] && begam < bg[Npt - 1]) interp = sp3->Eval(begam);
	delete sp3;
	return 100 * interp;
}
void Dndx(){
  dndx t;
  t.Loop();
}