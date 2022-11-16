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
   TFile *file=new TFile("dndx.root","recreate");
   TH1F *frac[25];       //each bin histogram
   for(Int_t i=0;i<25;i++){ char t_hname[25];
   sprintf(t_hname,"hg1%02d",i);
   frac[i]=new TH1F(t_hname,"",200,0,200);}
   TH1F *pth1[25];      
  for(Int_t i=0;i<25;i++){ char pthname[25];
    sprintf(pthname,"hpt1%02d",i);
    pth1[i]=new TH1F(pthname,"",100,pow(10,-0.2+i*0.09),pow(10,-0.2+(i+1)*0.09));}


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
         h->Fit("gaus","Q 0");//Q:Quiet mode (minimum printing)  0:拟合后不绘制直方图和拟合函数，但与选项"N"相反，它将拟合函数存储在函数直方图列表中
         TF1 *f =(TF1 *)h->GetFunction("gaus");
         sigma = f->GetParameter(2);
         // cout<<"sigma ="<<sigma<<endl;
         Double_t eff_w,sigma_w;
         eff_w = truth_ncl*0.01*(-0.007309)/(Track_L_DC[j]*TMath::Sqrt(1-pow(Track_Truth_CosTheta[j],2)))+1.245497;
         sigma_w = eff_w*TMath::Sqrt(truth_ncl);


         for ( int jp=0;jp<25;jp++){ if(bg>pow(10,-0.2+jp*0.09) && bg<pow(10,-0.2+(jp+1)*0.09))
         {frac[jp]->Fill(sigma);pth1[jp]->Fill(bg); } }/* } */

         // his_f->Fill(sigma/sigma_w);
         delete f;
         delete h;
         }
      }
   }
   
  double aver1[25];
  double pt1[25];
  for(int i=0;i<25;i++){
   aver1[i]=frac[i]->GetMean();
   frac[i]->Write();
   pt1[i]=pth1[i]->GetMean();
   cout<<"aver1 = "<<aver1[i]<<endl;
   cout<<"pth1 = "<<pt1[i]<<endl;
  }
   file->Close();

   TGraph *g = new TGraph(25,pt1,aver1);
   TCanvas *c1=new TCanvas("c1","c1",700,700);
   g->SetMarkerStyle(20);

   g->GetYaxis()->SetTitle("#sigma'");
   g->GetXaxis()->SetTitle("bg");
   g->Draw();
   c1->SetLogx();
   c1->Print("sigma_dndx.png");
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