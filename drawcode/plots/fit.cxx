#include<TH2.h>
#include<TStyle.h>
#include<TCanvas.h>
#include<cmath>
#include <iostream.h>
#include <iomanip.h>

void fit(){
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  //gStyle->SetTitleColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFillColor(0);

  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);
  
  gStyle->SetTitleFont(22,"xyz");  // set the all 3 axes title font
  gStyle->SetTitleFont(22," ");    // set the pad title font
  gStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  gStyle->SetLabelFont(22,"xyz");
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetTextFont(22);
  gStyle->SetTextSize(0.08);
  gStyle->SetStatFont(22);
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineWidth(2);  
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  
  gStyle->SetErrorX(0.001);
 
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gSystem->Load("libRooFit") ;
  using namespace RooFit ;

  float d_recoilmass,Likelihood;  

  TChain *chain1 = new TChain("TreeS","");
  chain1->Add("/scratchfs/bes/chenzx/cepc/TMVA/Output/hgg_BDT.root");
  chain1->SetBranchAddress("recoilmass", &d_recoilmass);
  chain1->SetBranchAddress("Likelihood", &Likelihood);

  TChain *chain2 = new TChain("TreeB","");
  chain2->Add("/scratchfs/bes/chenzx/cepc/TMVA/Output/hgg_BDT.root");
  chain2->SetBranchAddress("recoilmass", &d_recoilmass);
  chain2->SetBranchAddress("Likelihood", &Likelihood);

  RooRealVar x("x","M_{Higgs} GeV",120,150);
  //bw1 for etap
  RooRealVar* weight = new RooRealVar("weight", "weight", 0.0, 1000 );
  RooArgSet* ArgSet = new RooArgSet("args");
  ArgSet->add(x);
  ArgSet->add(*weight);
  Data = new RooDataSet("data", "A dataset", *ArgSet, "weight");

  Long64_t nevent1 = chain1->GetEntries();
  for(Long64_t j=0; j<nevent1; j++ )
  {
      chain1->GetEntry(j);

      if(!(d_recoilmass>120&&d_recoilmass<150&&Likelihood>0.11))continue;

      x.setVal(d_recoilmass);
      Data->add(*ArgSet, 1.0);
  }

  Long64_t nevent2 = chain2->GetEntries();
  for(Long64_t j=0; j<nevent2; j++ )
  {
      chain2->GetEntry(j);

      if(!(d_recoilmass>120&&d_recoilmass<150&&Likelihood>0.11))continue;

      x.setVal(d_recoilmass);
//      Data->add(*ArgSet, 1.22);
      Data->add(*ArgSet, 1.0);
  } 

//  RooRealVar mean("mean","mean",1.25227e+02,124.2,126);
  RooRealVar mean("mean","mean",1.25239e+02);
//  RooRealVar sigma("sigma","sigma",2.84170e-01,0.1,0.9);
  RooRealVar sigma("sigma","sigma",2.88009e-01);
//  RooRealVar a("a","coefficient a",-8.72979e-01,-10,0);
  RooRealVar a("a","coefficient a",-9.10730e-01);
//  RooRealVar n("n","coefficient n",1.27138e+00,0,10);
  RooRealVar n("n","coefficient n",1.04237e+00);
  RooCBShape sig1("sig1","sig1",x,mean,sigma,a,n);

//  RooRealVar a1("a1","coefficient a#1",-1.9449e-01,-10,10);
//  RooRealVar a2("a2","coefficient a#2",-9.1663e-02,-10,10);
//  RooRealVar a3("a3","coefficient a#3",-5.6504e-02,-10,10);


  RooRealVar a1("a1","coefficient a#1",-7.10465e-02);
  RooRealVar a2("a2","coefficient a#2",-5.39687e-02);
  RooRealVar a3("a3","coefficient a#3",-3.27216e-03);

  RooChebychev bkg("bkg","bkg",x,RooArgList(a1,a2,a3));  

  RooRealVar nsig1("nsig1","number of signal1 events",1.87378e+04,10000,50000);
  RooRealVar nbkg("nbkg","number of background events",13000,10000,300000);
//  RooRealVar nbkg("nbkg","number of background events",2.77319e+04);

  RooAddPdf model("model","sig1+bkg",RooArgList(sig1,bkg),RooArgList(nsig1,nbkg));
  //  RooAddPdf model("model","sig1",RooArgList(sig1),RooArgList(nsig1));

  RooPlot* xframe = x.frame(100);

  RooFitResult* RRR= model.fitTo(*data,Save(1));
  Double_t  aaa=RRR->minNll();
  
  data->plotOn(xframe);
  model.plotOn(xframe);
  cout<<"chisq of fit:"<<xframe->chiSquare(9)<<endl;
  cout << "-log(L) at minimum = "<<setprecision(8)<<aaa<< endl ;
  model.plotOn(xframe,Components(bkg),LineStyle(1),LineWidth(4),LineColor(3));
  model.plotOn(xframe,Components(sig1),LineStyle(1),LineWidth(4),LineColor(2));

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  xframe->Draw();
  xframe->GetXaxis()->SetTitle("M_{rec}(GeV/c^{2})");
  xframe->GetYaxis()->SetTitle("Events/(0.3GeV/c^{2})");
  xframe->GetXaxis()->SetLabelFont(22);
  xframe->GetYaxis()->SetLabelFont(22);
  xframe->GetXaxis()->CenterTitle();
  xframe->GetYaxis()->CenterTitle();
  xframe->GetXaxis()->SetNdivisions(404);
  xframe->GetYaxis()->SetNdivisions(404);
  xframe->GetYaxis()->SetRangeUser(0,8500);

  xframe->GetYaxis()->SetTitleOffset(1.10);
  xframe->GetYaxis()->SetLabelSize(0.05);
  xframe->GetYaxis()->SetLabelOffset(0.01);
  xframe->GetYaxis()->SetTitleSize(0.06);
  xframe->GetYaxis()->SetTitleFont(22);
  xframe->GetXaxis()->SetTitleOffset(1.10);
  xframe->GetXaxis()->SetLabelSize(0.05);
  xframe->GetXaxis()->SetLabelOffset(0.015);
  xframe->GetXaxis()->SetTitleSize(0.06);
  xframe->GetXaxis()->SetTitleFont(22);

  TLatex * CEPC = new TLatex(0.89,0.90, "CEPC");
  CEPC->SetNDC(); 
  CEPC->SetTextFont(22);
  CEPC->SetTextSize(0.08);
  CEPC->SetTextAlign(33);
  CEPC->Draw();
  
  TLatex * prelim = new TLatex(0.92,0.82, "Preliminary");
  prelim->SetNDC(); 
  prelim->SetTextFont(22);
  prelim->SetTextSize(0.055);
  prelim->SetTextAlign(33);
  prelim->Draw();
 
  TLine *xyz_1 = new TLine(130,7600,132,7600);
  xyz_1->SetLineColor(kBlue);
  xyz_1->SetLineWidth(4);
  xyz_1->Draw();
  
  TLine *xyz_2 = new TLine(130,7000,132,7000);
  xyz_2->SetLineColor(kRed);
  xyz_2->SetLineWidth(4);
  xyz_2->Draw();
  
  TLine *xyz_3 = new TLine(130,6400,132,6400);
  xyz_3->SetLineColor(kGreen);
  xyz_3->SetLineWidth(4);
  xyz_3->Draw();
  
 
  TLatex l;
  l.SetNDC();
  l.SetTextFont(22);
//  l.SetTextAlign(33);
  l.SetTextSize(0.035);
  l.DrawLatex(0.50,0.855,"Fit");
  l.DrawLatex(0.50,0.805,"Signal");
  l.DrawLatex(0.50,0.755,"Backgounds");

  c1->Print("figs/recoilfit.eps"); 
}
