#include <stdlib.h>
#include "TaJetProcessor.h"
#include <iostream>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TVector3.h"


// 090128 large change: angle based -> invmass based

using namespace lcio ;
using namespace marlin ;
using namespace std;

TaJetProcessor aTaJetProcessor ;

TaJetProcessor::TaJetProcessor() : Processor("TaJetProcessor")
{

	// modify processor description
	_description = "Tau Jet Finding" ;

	// register steering parameters: name, description, class-variable, default value
	registerInputCollection( LCIO::MCPARTICLE,
			"CollectionName" , 
			"Name of the MCParticle collection"  ,
			_colMCName ,
			std::string("MCParticle") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"PFOCollection" ,
			"PandoraPFA PFO collection",
			_pfoCollectionName,
			std::string("PandoraPFOs"));

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"OutputCollection",
			"TaJet output collection",
			_outputCollectionName,
			std::string("TaJets"));

	registerProcessorParameter( "MeasurementSmearing" ,
			"Additional acception angle for measurement smearing [rad]",
			_smearAngle,0.02);
	
	registerProcessorParameter( "LargestAngle",
			"Largest acceptance angle to be associated [rad]", _largestAngle,1.0);
	
	registerProcessorParameter( "TauMass", "Tau mass for jet finding [GeV] default=2.0", _tauMass,2.0);
	
	registerProcessorParameter( "MaxTracks", "Maximum number of tracks to be survive in preselection",_maxTracks,10);

}


void TaJetProcessor::init() { 

	streamlog_out(DEBUG) << "   init called  " 
		<< std::endl ;

	// usually a good idea to
	printParameters() ;

	_nRun = 0 ;
	_nEvt = 0 ;

}

void TaJetProcessor::processRunHeader( LCRunHeader* run) { 

	_nRun++ ;
} 

int compareEnergy(const void *pPart1, const void *pPart2)
{
	return int((*((ReconstructedParticle **)pPart2))->getEnergy()*10000 - (*((ReconstructedParticle **)pPart1))->getEnergy()*10000);
}

void TaJetProcessor::processEvent( LCEvent * evt ) { 

	// this gets called for every event 
	// usually the working horse ...


	// obtain Reconstructed Particles
	LCCollection* col = evt->getCollection(_pfoCollectionName);
	if(!col){streamlog_out(MESSAGE) << "Failed to obtain PFO collection!" << endl;return;}
	int nPart = col->getNumberOfElements();

	// output collection
	LCCollectionVec* pJetsCol= new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

	int ntrack = 0;
	for(int i=0; i < nPart; i++) {

		ReconstructedParticle *pPart = dynamic_cast<ReconstructedParticle*>( col->getElementAt( i ) );

		if(pPart->getCharge())ntrack ++;

		if(ntrack > _maxTracks){
			// skip event
			streamlog_out(DEBUG) << "Skip event: ntrack exceeds " << _maxTracks << std::endl;
			evt->addCollection(pJetsCol,_outputCollectionName);
			return;
		}
	}

	ReconstructedParticle **pPart = new ReconstructedParticle *[nPart];

	for(int i=0; i < nPart; i++) {
		pPart[i] = dynamic_cast<ReconstructedParticle*>( col->getElementAt( i ) );

		streamlog_out(MESSAGE) << "TAJPFO " << pPart[i]->getCharge() << " " << pPart[i]->getEnergy() << " " << pPart[i]->getMomentum()[2] << std::endl;
	}

	// sort by energy
	qsort(pPart,nPart,sizeof(ReconstructedParticle *),compareEnergy);

	streamlog_out(MESSAGE) << "TaJet reconstructed energies:" << endl;
	for(int i=0;i<nPart;i++){
		streamlog_out(MESSAGE) << pPart[i]->getEnergy() << (pPart[i]->getCharge() ? "C" : "N") << endl;
	}

	// Jet finding: around charged particle
	// because we do not care neutral jets for tau finding
	streamlog_out(MESSAGE) << "Tau reconstruction start." << endl;
	for(int i=0;i<nPart;i++)
	{
		ReconstructedParticle *pCur = pPart[i];
		if(pCur == NULL)continue;				// already associated particle skipped.
		if(!pCur->getCharge())continue;	// neutral particle skipped.

		streamlog_out(MESSAGE) << "Targetting to " << i << "th charged particle: energy = " << pCur->getEnergy() << endl;
		pPart[i] = NULL;	// detach from search list

		// create output jet
		ReconstructedParticleImpl *pJet = new ReconstructedParticleImpl;
		pJet->addParticle(pCur);

		double jetEnergy = pCur->getEnergy();
		TVector3 jetMomentum(pCur->getMomentum());
		double jetCharge = pCur->getCharge();

		streamlog_out(MESSAGE) << "Direction: ( " << pCur->getMomentum()[0] << ", " << pCur->getMomentum()[1] << ", "
			<< pCur->getMomentum()[2] << " )" << endl;


		int nPartAssoc = 1, nPartAssocPriv = 0;

		while(nPartAssoc != nPartAssocPriv){
			for(int j=0;j<nPart;j++)
			{
				ReconstructedParticle *pCur2 = pPart[j];
				if(pCur2 == NULL)continue;	 // already associated particle skipped.

				TVector3 vec2(pCur2->getMomentum());

				double angle = vec2.Angle(jetMomentum);

				// calculate invariant mass squared
				double invmass2 = pow(jetEnergy + pCur2->getEnergy(), 2) - (jetMomentum + vec2).Mag2();

				streamlog_out(MESSAGE) << "Candidate particle found: energy = " << pCur2->getEnergy() << ", " ;
				streamlog_out(MESSAGE) << "Direction: ( " << pCur2->getMomentum()[0] << ", " << pCur2->getMomentum()[1] << ", "
					<< pCur2->getMomentum()[2] << " ), charge = " << pCur2->getCharge()
					<< ", angle = " << angle << ", invmass = " << sqrt(invmass2) << endl;

				if(angle > _largestAngle) continue;	 // acceptance angle cut; skipped.
				if(invmass2 > _tauMass*_tauMass)continue;	// invmass cut; skipped.

				streamlog_out(MESSAGE) << "Candidate particle associated" << endl;
				pPart[j] = NULL;

				jetEnergy += pCur2->getEnergy();
				jetMomentum += vec2;
				jetCharge += pCur2->getCharge();
				pJet->addParticle(pCur2);

			}
			nPartAssocPriv = nPartAssoc;
			nPartAssoc = pJet->getParticles().size();
		}
		pJet->setEnergy(jetEnergy);
		double mom[3];
		mom[0] = jetMomentum.x();
		mom[1] = jetMomentum.y();
		mom[2] = jetMomentum.z();
		pJet->setMomentum(mom);
		pJet->setCharge(jetCharge);

		pJetsCol->addElement(pJet);

		streamlog_out(MESSAGE) << "TAJJETS " << jetEnergy << " " << jetMomentum.z() << std::endl;

	}
	// neutral particles
	for(int i=0;i<nPart;i++)
	{
		ReconstructedParticle *pCur = pPart[i];
		if(pCur == NULL)continue;				// already associated particle skipped.

		streamlog_out(MESSAGE) << "Remaining neutral particle " << i << ", " << pCur->getEnergy() << endl;
		pPart[i] = NULL;	// detach from search list

		// create output jet
		ReconstructedParticleImpl *pJet = new ReconstructedParticleImpl;
		pJet->addParticle(pCur);

		pJet->setEnergy(pCur->getEnergy());
		pJet->setMomentum(pCur->getMomentum());
		pJet->setCharge(pCur->getCharge());

		streamlog_out(MESSAGE) << "Direction: ( " << pCur->getMomentum()[0] << ", " << pCur->getMomentum()[1] << ", "
			<< pCur->getMomentum()[2] << " )" << endl;
		pJetsCol->addElement(pJet);
	}
	streamlog_out(MESSAGE) << "Tau reconstruction result: " << pJetsCol->getNumberOfElements() << " jets found." << endl;
	evt->addCollection(pJetsCol,_outputCollectionName);

	for(int i=0;i<pJetsCol->getNumberOfElements();i++)
	{
		ReconstructedParticle *ppart = (ReconstructedParticle *)pJetsCol->getElementAt(i);
		streamlog_out(MESSAGE) << "Jet " << i << ", charge: " << ppart->getCharge() << ", "
			<< "energy: " << ppart->getEnergy() << ", "
			<< "costheta: " << ppart->getMomentum()[2] / ppart->getEnergy() << std::endl;
	}

	//-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

	_nEvt ++ ;
}



void TaJetProcessor::check( LCEvent * evt ) { 
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TaJetProcessor::end(){ 

	//   std::cout << "TaJetProcessor::end()  " << name() 
	// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	// 	    << std::endl ;

}

