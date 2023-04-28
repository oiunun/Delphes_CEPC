#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>

#include <marlin/VerbosityLevels.h>
#include "ISRFinderProcessor.h"


using namespace lcio ;
using namespace marlin ;

ISRFinderProcessor aISRFinderProcessor ;

ISRFinderProcessor::ISRFinderProcessor()
	: Processor("ISRFinderProcessor") 
{

	// Processor description
	_description = "ISR Photon Finder Processor" ;

	// register steering parameters: name, description, class-variable, default value
	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputCollection" ,
			"Input collection of ReconstructedParticles",
			_inputPFOsCollection,
			std::string("PandoraPFOsWithoutIsoLep"));

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"OutputCollectionWithoutISR",
			"Copy of input collection but without the isolated leptons and initial state radiation",
			_outputPFOsRemovedISRCollection,
			std::string("PandoraPFOsWithoutISRWithoutIsoLep") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"OutputCollectionISR",
			"Output collection of initial state radiation",
			_outputISRCollection,
			std::string("ISR") );
}


void ISRFinderProcessor::init() { 
	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters() ;
}

void ISRFinderProcessor::processRunHeader( LCRunHeader* run) { 
} 

void ISRFinderProcessor::processEvent( LCEvent * evt ) { 

	_pfoCol = evt->getCollection( _inputPFOsCollection ) ;

	// Output PFOs removed isr 
	LCCollectionVec* otPFOsRemovedISRCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
	otPFOsRemovedISRCol->setSubset(true) ;

	// Output PFOs of isr
	LCCollectionVec* otISRCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	otISRCol->setSubset(true);

	// PFO loop
	int npfo = _pfoCol->getNumberOfElements();
	for (int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		if ( IsISR( pfo ) ) 
			otISRCol->addElement( pfo );
		else 
			otPFOsRemovedISRCol->addElement( pfo );
	}

	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

	// Add PFOs to new collection
	evt->addCollection( otPFOsRemovedISRCol, _outputPFOsRemovedISRCollection.c_str() );
	evt->addCollection( otISRCol, _outputISRCollection.c_str() );
}

void ISRFinderProcessor::check( LCEvent * evt ) { 
}

void ISRFinderProcessor::end() { 
}

bool ISRFinderProcessor::IsCharged( ReconstructedParticle* pfo ) {
	if ( pfo->getCharge() == 0 ) return true;
	return false;
}

bool ISRFinderProcessor::IsEnergyandCos( ReconstructedParticle* pfo ){
	TVector3 moment(pfo->getMomentum());
	double cos = moment.CosTheta();
	if ( (pfo->getEnergy() > 80 && fabs(cos) > 0.85) || (pfo->getEnergy() > 30 && fabs(cos) > 0.95)) return true;
	return false;
}

bool ISRFinderProcessor::PhotonSelection( ReconstructedParticle* pfo) {
	float ecal = 0;
	float hcal = 0;
	float cale[2];
	double ratio = 0;
	std::vector<lcio::Cluster*> clusters = pfo->getClusters();
	for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
			iCluster!=clusters.end();
			++iCluster) {
		ecal += (*iCluster)->getSubdetectorEnergies()[0];
		hcal += (*iCluster)->getSubdetectorEnergies()[1];
	}
	cale[0] = ecal;
	cale[1] = hcal;
	ratio = ecal/(hcal+ecal);
	if (ratio > 0.9) return true;
	return false;
}

bool ISRFinderProcessor::IsISR( ReconstructedParticle* pfo ) {
	if ( !IsCharged(pfo) )
		return false;

	if ( !PhotonSelection(pfo) )
		return false;
	
	if ( !IsEnergyandCos(pfo) )
		return false;

	return true;
}

