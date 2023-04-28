// *****************************************************
// Processor for overlay optimization
//                        ----Junping
// *****************************************************
#include "overlayOptimizationProcessor.h"
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


overlayOptimizationProcessor aoverlayOptimizationProcessor ;


overlayOptimizationProcessor::overlayOptimizationProcessor() : Processor("overlayOptimizationProcessor") {

	// modify processor description
	_description = "overlayOptimizationProcessor does whatever it does ..." ;


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
			"OutputPFOsAfterFiltCollection",
			"Name of the PFOs collection selected angle cut",
			_colNewPFOs,
			std::string("FiltedPFOs") );


	registerProcessorParameter("CosThetaCut",   
			"remove small angle PFOs",  
			_CosThetaCut,   
			double(1.00)  );


}

void overlayOptimizationProcessor::init() { 

	streamlog_out(DEBUG) << "   init called  " 
		<< std::endl ;


	// usually a good idea to
	printParameters() ;

	_nRun = 0 ;
	_nEvt = 0 ;


}

void overlayOptimizationProcessor::processRunHeader( LCRunHeader* run) { 

	_nRun++ ;
} 

void overlayOptimizationProcessor::processEvent( LCEvent * evt ) { 


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
		//Double_t energy = recPart->getEnergy();
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
		TVector3 momentum = TVector3(recPart->getMomentum());
		Double_t cosTheta = momentum.CosTheta();
		if ( fabs(cosTheta) < _CosThetaCut ) newPFOs.push_back(recPart);
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

void overlayOptimizationProcessor::check( LCEvent * evt ) { 
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}

void overlayOptimizationProcessor::end(){ 

	cerr << "overlayOptimizationProcessor::end()  " << name() 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< endl ;
}
