/*
root -l examples/Example1.C'("delphes_output.root")'
*/

//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "../classes/DelphesClasses.h"
#include "../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
//#endif
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

#include "cepcplotstyle.C"


using namespace RooFit;
using namespace std   ;

//------------------------------------------------------------------------------

void FitRecoil_mu2(const char * inputFile )
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  float d_recoilmass;
  Int_t number_select = 0;  

  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TLorentzVector cms(0,0,0,240);
  // Book histograms
  TH1 *histMass = new TH1F("mass", "M_{inv}(#mu_{1},#mu_{2})", 100, 75, 110);


	RooRealVar x("x","M_{Higgs} GeV",120,140);
	RooRealVar* weight = new RooRealVar("weight", "weight", 0.0, 1000 );
	RooArgSet* ArgSet = new RooArgSet("args");
	ArgSet->add(x);
	ArgSet->add(*weight);
	RooAbsData* data = new RooDataSet("data", "A dataset", *ArgSet);

  // Loop over all events
  for(Int_t entry = 0; entry < 100000; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    Muon *mu1 ,*mu2;
    // If event contains at least 2 electrons
    if(branchMuon->GetEntries() > 1){  
      // Take correct two muons
      mu1 = (Muon *) branchMuon->At(0);
      mu2 = (Muon *) branchMuon->At(1);
      number_select ++;
      }
      // Plot their invariant mass
      d_recoilmass = ((cms - mu1->P4() - mu2->P4())).M();
      histMass->Fill(d_recoilmass);
 		if(!(d_recoilmass>120&&d_recoilmass<140))continue;
		x.setVal(d_recoilmass);
		data->add(*ArgSet,1.0);
     
    }

	//use a function which is composed of gaus and exponential 
	//function to fit the two muons' invariant mass spectrum
	
	RooRealVar mean("mean","mean",1.25239e+02,124,126);
	RooRealVar sigma("sigma","sigma",2.88009e-01,0.1,0.25);
	RooRealVar a("a","coefficient a",-9.10730e-01);
	RooRealVar n("n","coefficient n",1.04237e+00);
	RooCBShape sig1("sig1","Signal",x,mean,sigma,a,n);

	RooRealVar nsig1("nsig1","number of signal1 events",1.60378e+04,0,100000);
	RooAddPdf  model("model","S Fit",RooArgList(sig1), RooArgList(nsig1));

	vector<RooAbsPdf*> pdfList; pdfList.clear(); 
	RooFitResult* RRR= model.fitTo(*data, Verbose(kFALSE), Extended(kTRUE), Save(1), PrintLevel(1));
	RRR->Print("v");

	pdfList.push_back(&model);

	PlotDataRooFit("fig/mass_recoil", data, x, pdfList, "M_{recoil}^{#mu^{+}#mu^{-}}[GeV]", "Entries/0.2 GeV", 
			true,
			"Z#rightarrow #mu^{+}#mu^{-}",
			-1 // maximumm of y-axis, defaut -1
			);

  cout<<"numberOfEntries = "<<numberOfEntries<<endl;
  cout<<"number_select = "<<number_select<<endl;
  

}

