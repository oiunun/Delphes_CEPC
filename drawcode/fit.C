#include <cmath>
#include <iostream>
#include <iomanip>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <RooChebychev.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>

#include <cepcplotstyle.C>


using namespace RooFit;
using namespace std   ;

void fit(){

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
	RooRealVar* weight = new RooRealVar("weight", "weight", 0.0, 1000 );
	RooArgSet* ArgSet = new RooArgSet("args");
	ArgSet->add(x);
	ArgSet->add(*weight);
	RooAbsData* data = new RooDataSet("data", "A dataset", *ArgSet, "weight");

	Long64_t nevent1 = chain1->GetEntries();
	for(Long64_t j=0; j<nevent1; j++ )
	{
		chain1->GetEntry(j);
		if(!(d_recoilmass>120&&d_recoilmass<150&&Likelihood>0.11))continue;
		x.setVal(d_recoilmass);
		data->add(*ArgSet, 1.0);
	}

	Long64_t nevent2 = chain2->GetEntries();
	for(Long64_t j=0; j<nevent2; j++ )
	{
		chain2->GetEntry(j);
		if(!(d_recoilmass>120&&d_recoilmass<150&&Likelihood>0.11))continue;
		x.setVal(d_recoilmass);
		data->add(*ArgSet, 1.0);
	} 

	RooRealVar mean("mean","mean",1.25239e+02);
	RooRealVar sigma("sigma","sigma",2.88009e-01);
	RooRealVar a("a","coefficient a",-9.10730e-01);
	RooRealVar n("n","coefficient n",1.04237e+00);
	RooCBShape sig1("sig1","Signal",x,mean,sigma,a,n);

	RooRealVar a1("a1","coefficient a#1",-7.10465e-02);
	RooRealVar a2("a2","coefficient a#2",-5.39687e-02);
	RooRealVar a3("a3","coefficient a#3",-3.27216e-03);
	RooChebychev bkg("bkg","Background",x,RooArgList(a1,a2,a3));  

	RooRealVar nsig1("nsig1","number of signal1 events",1.87378e+04,10000,50000);
	RooRealVar  nbkg("nbkg" ,"number of background events",13000,10000,300000);
	RooAddPdf  model("model","S+B Fit", RooArgList(sig1,bkg), RooArgList(nsig1,nbkg));

	vector<RooAbsPdf*> pdfList; pdfList.clear(); 
	RooFitResult* RRR= model.fitTo(*data, Verbose(kFALSE), Extended(kTRUE), Save(1), PrintLevel(1));
	RRR->Print("v");

	pdfList.push_back(&model);
	pdfList.push_back(&sig1 );
	pdfList.push_back(&bkg  );

	PlotDataRooFit("fit/recoil", data, x, pdfList, "M_{recoil}^{#mu^{+}#mu^{-}}[GeV]", "Entries/0.3 GeV", 
			true,
			"Z#rightarrow #mu^{+}#mu^{-}, H#rightarrow b#bar{b}; #int Ldt = 5 ab^{-1}",
			5500 // maximumm of y-axis, defaut -1
			);
}
