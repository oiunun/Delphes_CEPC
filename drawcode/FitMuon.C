/*
root -l examples/Example1.C'("delphes_output.root")'
*/

//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "../classes/DelphesClasses.h"
#include "../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
//#endif
#include "TMath.h"
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
#include <RooKeysPdf.h>
#include <RooFitResult.h>
#include <RooFFTConvPdf.h>

#include "cepcplotstyle.C"


using namespace RooFit;
using namespace std   ;

//------------------------------------------------------------------------------

void FitMuon(const char * inputFile )
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  float d_mass;
  Int_t number_select = 0;  

  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  // Book histograms

	RooRealVar x("x","M_{Z} GeV",75,105);

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
      
      // Plot their invariant mass
    d_mass = ((  mu1->P4() +  mu2->P4())).M();
 		if(!(d_mass>75&&d_mass<105))continue;
		x.setVal(d_mass);
		data->add(*ArgSet,1.0);
     }
    }

	//use a function which is composed of gaus and exponential 
	//function to fit the two muons' invariant mass spectrum

	RooRealVar GaussM("GaussM","m",   91  , 89, 93 );
	RooRealVar GaussS("GaussS","c10",   0.6 ,  0, 5);//,-5e-3,5e-3);
	// RooGaussian Gauss("Gauss","gauss", x, GaussM, GaussS);


	RooRealVar M2("M2","m2",  91.1876 , 90.8 ,91.5 );
	RooRealVar S2("S2","c102", 2 , 0,2.5);//,-5e-3,5e-3);


  // RooKeysPdf sig1x("sig1x", "sig1x", x, data, RooKeysPdf::NoMirror,1);
  // RooGenericPdf sig1x("sigx1","1.2476/(3.1415926*(pow((@0-91.18),2)+1.5565))",RooArgList(x));
  // RooBreitWigner sig1x("sig1x","sig1x",x,M2,S2);

  // RooNumConvPdf sig("sig","",x,sig1x,Gauss);
  RooVoigtian sig("sig","",x,M2,S2,GaussS);
	RooRealVar nsig1("nsig1","number of signal1 events",17041,0,200000);

	RooAddPdf  model("model"," Fit",RooArgList(sig), RooArgList(nsig1));

  // TCanvas *c1=new TCanvas("c1","c1",1300,700);
  // sig.fitTo(*data);
  // RooPlot *frame1 = x.frame(Title("Cyclical convolution in angle psi"));
  // data->plotOn(frame1);
  // sig.plotOn(frame1);
  // c1->Print("../fig/mass_muon.png");
  
	vector<RooAbsPdf*> pdfList; pdfList.clear(); 
  x.setRange("R",89,105);
	RooFitResult* RRR= model.fitTo(*data, Verbose(kFALSE), Range("R"), Extended(kTRUE), Save(1), PrintLevel(1));
	RRR->Print("v");

	pdfList.push_back(&model);

	PlotDataRooFit("fig/mass_muon", data, x, pdfList, "M^{#mu^{+}#mu^{-}}[GeV]", "Entries/0.3 GeV", 
			true,
			"Z#rightarrow #mu^{+}#mu^{-}",
			-1);

  cout<<"numberOfEntries = "<<numberOfEntries<<endl;
  cout<<"number_select = "<<number_select<<endl;
  

}

