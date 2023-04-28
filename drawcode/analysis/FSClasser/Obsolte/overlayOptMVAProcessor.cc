// *****************************************************
// Processor for overlay optimization
//                        ----Junping
// *****************************************************
#include "overlayOptMVAProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TMVA/Reader.h"

#include "Utilities.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace TMVA;

overlayOptMVAProcessor aoverlayOptMVAProcessor ;


overlayOptMVAProcessor::overlayOptMVAProcessor() : Processor("overlayOptMVAProcessor") {

	// modify processor description
	_description = "overlayOptMVAProcessor does whatever it does ..." ;


	// register steering parameters: name, description, class-variable, default value

	registerInputCollection( LCIO::MCPARTICLE,
			"InputMCParticlesCollection" , 
			"Name of the MCParticle collection"  ,
			_colMCP ,
			std::string("MCParticlesSkimmed") ) ;

	registerInputCollection( LCIO::LCRELATION,
			"InputMCTruthLinkCollection" , 
			"Name of the MCTruthLink collection"  ,
			_colMCTL ,
			std::string("RecoMCTruthLink") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputPandoraPFOsCollection" , 
			"Name of the PandoraPFOs collection"  ,
			_colPFOs ,
			std::string("PandoraPFOs") ) ;

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputPFOsAfterMVACollection",
			"Name of the PFOs collection selected by MVA",
			_colNewPFOs,
			std::string("newpfos") );

	registerProcessorParameter("DirOfOverlayWeights1",
			"Directory of Weights for the overlay MVA Classification"  ,
			_overlay_weights1 ,
			std::string("overlay1_hww_weights") ) ;

	registerProcessorParameter("DirOfOverlayWeights2",
			"Directory of Weights for the overlay MVA Classification"  ,
			_overlay_weights2 ,
			std::string("overlay2_hww_weights") ) ;

	registerProcessorParameter("CutOnTheOverlayMVAOUTPUT1",
			"Cut on the 1st mva output"  ,
			_overlay_mva_cut1 ,
			float(0.0) ) ;

	registerProcessorParameter("CutOnTheOverlayMVAOUTPUT2",
			"Cut on the 2nd mva output"  ,
			_overlay_mva_cut2 ,
			float(0.0) ) ;

	registerProcessorParameter("CosThetaCut",   
			"remove small angle PFOs",  
			_CosThetaCut,   
			double(1.00)  );


}

void overlayOptMVAProcessor::init() { 

	streamlog_out(DEBUG) << "   init called  " 
		<< std::endl ;


	// usually a good idea to
	printParameters() ;

	_nRun = 0 ;
	_nEvt = 0 ;

	//  TFile *outRootFile = new TFile("output.root","RECREATE");

	// add variables
	for (Int_t i=0;i<2;i++) {
		TMVA::Reader *reader = new TMVA::Reader( "Color:Silent" );    
		reader->AddVariable( "pt",          &_pt );
		reader->AddVariable( "rapidity",    &_rapidity );
		if (i == 1) {
			reader->AddVariable( "z0",       &_z0 );
		}
		//      reader->AddVariable( "costheta",    &_costheta );
		// book the reader (method, weights)
		TString dir    = _overlay_weights1;
		if (i == 1) {
			dir = _overlay_weights2;
		}
		TString prefix = "TMVAClassification";
		TString methodName = "BDT method";
		TString weightfile = dir + "/" + prefix + "_" + "BDT.weights.xml";
		reader->BookMVA( methodName, weightfile ); 
		_readers.push_back(reader);
	}

}

void overlayOptMVAProcessor::processRunHeader( LCRunHeader* run) { 

	_nRun++ ;
} 

void overlayOptMVAProcessor::processEvent( LCEvent * evt ) { 


	// this gets called for every event 
	// usually the working horse ...
	_nEvt++;

	//cerr << endl << "Hello, Overlay MVA Selection!" << endl;

	// -- Read out PFO information --
	LCCollection *colPFO = evt->getCollection(_colPFOs);
	if (!colPFO) {
		std::cerr << "No PFO Collection Found!" << std::endl;
		throw marlin::SkipEventException(this);
	}

	Int_t nPFOs = colPFO->getNumberOfElements();
	//  cerr << "Number of PFOs: " << nPFOs << endl;
	LCCollectionVec *pNewPFOsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	pNewPFOsCollection->setSubset(true);
	std::vector<lcio::ReconstructedParticle*> newPFOs;

	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		Double_t energy = recPart->getEnergy();
		//Double_t charge = recPart->getCharge();
		TrackVec tckvec = recPart->getTracks();
		Int_t ntracks = tckvec.size();
		Double_t d0=0.,z0=0.,deltad0=0.,deltaz0=0.,nsigd0=0.,nsigz0=0.;
		if (ntracks > 0) {
			d0 = tckvec[0]->getD0();
			z0 = tckvec[0]->getZ0();
			deltad0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
			deltaz0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
			nsigd0 = d0/deltad0;
			nsigz0 = z0/deltaz0;
		}
		//Double_t r0 = TMath::Sqrt(d0*d0+z0*z0);
		//Double_t nsigr0 = TMath::Sqrt(nsigd0*nsigd0+nsigz0*nsigz0);
		TVector3 momentum = TVector3(recPart->getMomentum());
		//Double_t momentumMagnitude = momentum.Mag();
		Double_t pT = momentum.Pt();
		Double_t pZ = momentum.Pz();
		Double_t cosTheta = momentum.CosTheta();
		Double_t rapidity = TMath::Log((energy+pZ)/(energy-pZ))/2.;
		// evaluate the overlay neural-net output for each pfo
		_pt = pT;
		_rapidity = rapidity;
		_costheta = cosTheta;
		_z0 = z0;
		Double_t mva_overlay1 = -1., mva_overlay2 = -1.;
		mva_overlay1 = _readers[0]->EvaluateMVA( "BDT method"   );
		mva_overlay2 = _readers[1]->EvaluateMVA( "BDT method"   );
		//if(((charge==0||TMath::Abs(z0) > 1.) && mva_overlay1 > _overlay_mva_cut1)||(charge!=0&&TMath::Abs(z0)<1.&&mva_overlay2>_overlay_mva_cut2))
		if ( 1 )
		{
			if(fabs(cosTheta)<_CosThetaCut ) newPFOs.push_back(recPart);
		}
	}

	// save the other PFOs to a new collection
	for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=newPFOs.begin();iObj<newPFOs.end();++iObj) {
		pNewPFOsCollection->addElement(*iObj);
	}
	evt->addCollection(pNewPFOsCollection,_colNewPFOs.c_str());

	//-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

}



void overlayOptMVAProcessor::check( LCEvent * evt ) { 
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void overlayOptMVAProcessor::end(){ 

	cerr << "overlayOptMVAProcessor::end()  " << name() 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< endl ;

}
