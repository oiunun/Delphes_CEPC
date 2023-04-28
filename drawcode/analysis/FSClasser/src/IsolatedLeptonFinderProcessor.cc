#include "IsolatedLeptonFinderProcessor.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace lcio ;
using namespace marlin ;

LGISOlatedLeptonFinderProcessor aISOlatedLeptonFinderProcessor ;

LGISOlatedLeptonFinderProcessor::LGISOlatedLeptonFinderProcessor()
	: Processor("LGISOlatedLeptonFinderProcessor") {

		// Processor description
		_description = "ISOlated Lepton Finder Processor" ;

		// register steering parameters: name, description, class-variable, default value
		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"InputCollection" ,
				"Input collection of ReconstructedParticles",
				_inputPFOsCollection,
				std::string("ArborPFOs"));

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionWithoutIsolatedLepton",
				"Copy of input collection but without the isolated leptons",
				_outputPFOsRemovedIsoLepCollection,
				std::string("ArborPFOsWithoutIsoLep") );

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionIsolatedLeptons",
				"Output collection of isolated leptons",
				_outputIsoLepCollection,
				std::string("Isolep") );
		
		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionWithoutIsolatedZLepton",
				"Copy of input collection but without the isolated leptons from Z",
				_outputPFOsRemovedIsoZLepCollection,
				std::string("ArborPFOsWithoutIsoZLep") );

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollectionIsolatedZLeptons",
				"Output collection of isolated leptons from Z",
				_outputIsoZLepCollection,
				std::string("IsoZlep") );

		registerProcessorParameter( "CosConeAngle",
				"Cosine of the half-angle of the cone used in isolation criteria",
				_cosConeAngle,
				double(0.98));

		registerProcessorParameter( "UsePID",
				"Use primitive particle ID based on calorimeter energy deposits",
				_usePID,
				bool(true));

		registerProcessorParameter( "ElectronMinEnergyDepositByMomentum",
				"Electron ID: Minimum energy deposit divided by momentum",
				_electronMinEnergyDepositByMomentum,
				double(0.7));

		registerProcessorParameter( "ElectronMaxEnergyDepositByMomentum",
				"Electron ID: Maximum energy deposit divided by momentum",
				_electronMaxEnergyDepositByMomentum,
				double(1.4));

		registerProcessorParameter( "ElectronMinEcalToHcalFraction",
				"Electron ID: Minimum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_electronMinEcalToHcalFraction,
				double(0.9));

		registerProcessorParameter( "ElectronMaxEcalToHcalFraction",
				"Electron ID: Maximum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_electronMaxEcalToHcalFraction,
				double(1.0));

		registerProcessorParameter( "MuonMinEnergyDepositByMomentum",
				"Muon ID: Minimum energy deposit divided by momentum",
				_muonMinEnergyDepositByMomentum,
				double(0.0));

		registerProcessorParameter( "MuonMaxEnergyDepositByMomentum",
				"Muon ID: Maximum energy deposit divided by momentum",
				_muonMaxEnergyDepositByMomentum,
				double(0.3));

		registerProcessorParameter( "MuonMinEcalToHcalFraction",
				"Muon ID: Minimum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_muonMinEcalToHcalFraction,
				double(0.0));

		registerProcessorParameter( "MuonMaxEcalToHcalFraction",
				"Muon ID: Maximum Ecal deposit divided by sum of Ecal and Hcal deposits",
				_muonMaxEcalToHcalFraction,
				double(0.4));

		registerProcessorParameter( "UseImpactParameter",
				"Use impact parameter cuts for consistency with primary/secondary track",
				_useImpactParameter,
				bool(true));

		registerProcessorParameter( "ImpactParameterMinD0",
				"Minimum d0 impact parameter",
				_minD0,
				double(0.0));

		registerProcessorParameter( "ImpactParameterMaxD0",
				"Maximum d0 impact parameter",
				_maxD0,
				double(1e20));

		registerProcessorParameter( "ImpactParameterMinZ0",
				"Minimum z0 impact parameter",
				_minZ0,
				double(0.0));

		registerProcessorParameter( "ImpactParameterMaxZ0",
				"Maximum z0 impact parameter",
				_maxZ0,
				double(1e20));

		registerProcessorParameter( "ImpactParameterMin3D",
				"Minimum impact parameter in 3D",
				_minR0,
				double(0.0));

		registerProcessorParameter( "ImpactParameterMax3D",
				"Maximum impact parameter in 3D",
				_maxR0,
				double(0.01));

		registerProcessorParameter( "UseImpactParameterSignificance",
				"Use impact parameter significance cuts for consistency with primary/secondary track",
				_useImpactParameterSignificance,
				bool(false));

		registerProcessorParameter( "ImpactParameterMinD0Significance",
				"Minimum d0 impact parameter significance",
				_minD0Sig,
				double(0.0));

		registerProcessorParameter( "ImpactParameterMaxD0Significance",
				"Maximum d0 impact parameter significance",
				_maxD0Sig,
				double(1e20));

		registerProcessorParameter( "ImpactParameterMinZ0Significance",
				"Minimum z0 impact parameter significance",
				_minZ0Sig,
				double(0.0));

		registerProcessorParameter( "ImpactParameterMaxZ0Significance",
				"Maximum z0 impact parameter significance",
				_maxZ0Sig,
				double(1e20));

		registerProcessorParameter( "ImpactParameterMin3DSignificance",
				"Minimum impact parameter significance in 3D",
				_minR0Sig,
				double(0.0));

		registerProcessorParameter( "ImpactParameterMax3DSignificance",
				"Maximum impact parameter significance in 3D",
				_maxR0Sig,
				double(1e20));

		registerProcessorParameter( "UseRectangularIsolation",
				"Use rectangular cuts on track and cone energy",
				_useRectangularIsolation,
				bool(true));

		registerProcessorParameter( "IsolationMinimumTrackEnergy",
				"Minimum track energy for isolation requirement",
				_isoMinTrackEnergy,
				double(2));

		registerProcessorParameter( "ZIsolationMinimumTrackEnergy",
				"Minimum track energy for isolation requirement",
				_ZisoMinTrackEnergy,
				double(15));

		registerProcessorParameter( "IsolationMaximumTrackEnergy",
				"Maximum track energy for isolation requirement",
				_isoMaxTrackEnergy,
				double(1e20));

		registerProcessorParameter( "IsolationMinimumConeEnergy",
				"Minimum cone energy for isolation requirement",
				_isoMinConeEnergy,
				double(0));

		registerProcessorParameter( "IsolationMaximumConeEnergy",
				"Maximum cone energy for isolation requirement",
				_isoMaxConeEnergy,
				double(1e20));

		registerProcessorParameter( "UsePolynomialIsolation",
				"Use polynomial cuts on track and cone energy",
				_usePolynomialIsolation,
				bool(true));

		registerProcessorParameter( "IsolationPolynomialCutA",
				"Polynomial cut (A) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
				_isoPolynomialA,
				double(0));

		registerProcessorParameter( "IsolationPolynomialCutB",
				"Polynomial cut (B) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
				_isoPolynomialB,
				double(20));

		registerProcessorParameter( "IsolationPolynomialCutC",
				"Polynomial cut (C) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
				_isoPolynomialC,
				double(-300));

		registerProcessorParameter( "UseJetIsolation",
				"Use jet-based isolation",
				_useJetIsolation,
				bool(false));

		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollection" ,
				"Input collection of jets for isolation",
				_jetCollectionName,
				std::string("JetsForIsolation"));

		registerProcessorParameter( "JetIsolationVetoMinimumXt",
				"Minimum Xt in jet-based isolation",
				_jetIsoVetoMinXt,
				double(0.));

		registerProcessorParameter( "JetIsolationVetoMaximumXt",
				"Maximum Xt in jet-based isolation",
				_jetIsoVetoMaxXt,
				double(0.25));

		registerProcessorParameter( "JetIsolationVetoMinimumZ",
				"Mininum Z in jet-based isolation",
				_jetIsoVetoMinZ,
				double(0.));

		registerProcessorParameter( "JetIsolationVetoMaximumZ",
				"Maximum Z in jet-based isolation",
				_jetIsoVetoMaxZ,
				double(0.6));

		registerProcessorParameter( "mBoson",
				"Higgs or Z mass with lepton pair ",
				_mBoson,
				double(91.18));
	}


void LGISOlatedLeptonFinderProcessor::init() { 
	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters() ;
}

void LGISOlatedLeptonFinderProcessor::processRunHeader( LCRunHeader* run) { 
} 

void LGISOlatedLeptonFinderProcessor::processEvent( LCEvent * evt ) { 

	_rpJetMap.clear();

	_pfoCol = evt->getCollection( _inputPFOsCollection ) ;

	// Output PFOs removed isolated leptons 
	LCCollectionVec* otPFOsRemovedIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
	otPFOsRemovedIsoLepCol->setSubset(true) ;
	LCCollectionVec* otPFOsRemovedIsoZLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
	otPFOsRemovedIsoZLepCol->setSubset(true) ;

	// Output PFOs of isolated leptons
	LCCollectionVec* otIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	otIsoLepCol->setSubset(true);
	LCCollectionVec* otIsoZLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	otIsoZLepCol->setSubset(true);

	// Prepare jet/recoparticle map for jet-based isolation
	if (_useJetIsolation) {
		LCCollection *colJet = evt->getCollection(_jetCollectionName);
		int njet = colJet->getNumberOfElements();
		for (int i=0; i<njet; ++i) {
			ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( colJet->getElementAt(i) );
			for (ReconstructedParticleVec::const_iterator iter = jet->getParticles().begin();
					iter != jet->getParticles().end(); ++iter) {
				_rpJetMap.insert( std::make_pair( *iter, jet ) );
			}
		}
	}

	// PFO loop
	int npfo = _pfoCol->getNumberOfElements();
	for (int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		if ( IsISOlatedLepton( pfo ) ) 
			otIsoLepCol->addElement( pfo );
		else 
			otPFOsRemovedIsoLepCol->addElement( pfo );
	}

	// Lepton loop
	int nlep = otIsoLepCol->getNumberOfElements();
	vector<ReconstructedParticle*> ZLeps; 
	if( nlep>=2){
		
		ZLeps = getZLeptonPair( otIsoLepCol ) ;
	
		for (int i = 0; i < npfo; i++ ) {
			ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

			vector<ReconstructedParticle*>::iterator PosLep = find( ZLeps.begin(), ZLeps.end(), pfo );
			if ( ZLeps.end() != PosLep ) 
				otIsoZLepCol->addElement( pfo );
			else 
				otPFOsRemovedIsoZLepCol->addElement( pfo );
		}	
	
	}

	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

	// Add PFOs to new collection
	evt->addCollection( otPFOsRemovedIsoLepCol,  _outputPFOsRemovedIsoLepCollection.c_str() );
	evt->addCollection( otIsoLepCol,             _outputIsoLepCollection.c_str() );
	evt->addCollection( otPFOsRemovedIsoZLepCol, _outputPFOsRemovedIsoZLepCollection.c_str() );
	evt->addCollection( otIsoZLepCol,            _outputIsoZLepCollection.c_str() );
}

void LGISOlatedLeptonFinderProcessor::check( LCEvent * evt ) { 
}

void LGISOlatedLeptonFinderProcessor::end() { 
}

vector<ReconstructedParticle*> LGISOlatedLeptonFinderProcessor::getZLeptonPair( LCCollection* leps ){
	int nlep = leps->getNumberOfElements();
	double dMass=9999.0;
	vector<ReconstructedParticle*> ret; 
	for ( int i=0; i<nlep-1 ; i++){
		ReconstructedParticle* enflow1 = dynamic_cast<ReconstructedParticle*>(leps->getElementAt( i ));
		int pdgid1=enflow1->getType();
		Double_t energy1   = enflow1->getEnergy();
		TVector3 momentum1 = TVector3(enflow1->getMomentum());
		TLorentzVector p41 = TLorentzVector(momentum1,energy1);
		if( energy1<_ZisoMinTrackEnergy) continue;
		//
		for ( int j=i+1; j<nlep ; j++){
			ReconstructedParticle* enflow2 = dynamic_cast<ReconstructedParticle*>(leps->getElementAt( j ));
			int pdgid2=enflow2->getType();
			if(pdgid1+pdgid2!=0) continue;
			Double_t energy2   = enflow2->getEnergy();
			TVector3 momentum2 = TVector3(enflow2->getMomentum());
			TLorentzVector p42 = TLorentzVector(momentum2,energy2);
			if( energy2<_ZisoMinTrackEnergy) continue;
			//
			double dmass = fabs( (p41+p42).M() - _mBoson ); 
			if ( dmass < dMass ){
				dMass=dmass;
				ret.clear();
				ret.push_back(enflow1);
				ret.push_back(enflow2);
			}
		}
	}
	return ret;
}

bool LGISOlatedLeptonFinderProcessor::IsCharged( ReconstructedParticle* pfo ) {
	if ( pfo->getCharge() == 0 ) return false;
	return true;
}

bool LGISOlatedLeptonFinderProcessor::IsLepton( ReconstructedParticle* pfo ) {
	/*   
	double CalE[2];
	getCalEnergy( pfo , CalE );
	double ecale  = CalE[0];
	double hcale  = CalE[1];
	double p      = TVector3( pfo->getMomentum() ).Mag();
	double calByP = p>0 ? (ecale + hcale)/p : 0;
	double calSum = ecale+hcale;
	double ecalFrac = calSum>0 ? ecale / calSum : 0;

	//printf("ecal, hcal, p and pid = %8.3f %8.3f %8.3f %8d\n", ecale, hcale, p, pfo->getType());
	// electron
	if ( calByP >= _electronMinEnergyDepositByMomentum
			&& calByP <= _electronMaxEnergyDepositByMomentum
			&& ecalFrac >= _electronMinEcalToHcalFraction
			&& ecalFrac <= _electronMaxEcalToHcalFraction ){
		return true;
	}

	// muon
	if ( calByP >= _muonMinEnergyDepositByMomentum
			&& calByP <= _muonMaxEnergyDepositByMomentum
			&& ecalFrac >= _muonMinEcalToHcalFraction
			&& ecalFrac <= _muonMaxEcalToHcalFraction ){
		return true;
	}
	cout<<" type = " << pfo->getType()<<endl;
	*/
	if( abs(pfo->getType())==11 || abs(pfo->getType())==13 ) return true;

	return false;
}

bool LGISOlatedLeptonFinderProcessor::IsISOlatedLepton( ReconstructedParticle* pfo ) {

	if ( !IsCharged(pfo) )
		return false;

	if ( _usePID && !IsLepton(pfo) )
		return false;

	if ( _useImpactParameter && !PassesImpactParameterCuts(pfo) )
		return false ;

	if ( _useImpactParameterSignificance && !PassesImpactParameterSignificanceCuts(pfo) )
		return false ;

	if ( _useRectangularIsolation && !IsISOlatedRectangular(pfo) )
		return false;

	if ( _usePolynomialIsolation && !IsISOlatedPolynomial(pfo) )
		return false;

	if ( _useJetIsolation && !IsISOlatedJet(pfo) )
		return false;

	return true;
}

bool LGISOlatedLeptonFinderProcessor::IsISOlatedRectangular( ReconstructedParticle* pfo ) {

	double E     = pfo->getEnergy() ;
	double coneE = getConeEnergy( pfo )/E;

	if (E < _isoMinTrackEnergy) return false;
	if (E > _isoMaxTrackEnergy) return false;
	if (coneE < _isoMinConeEnergy) return false;
	if (coneE > _isoMaxConeEnergy) return false;

	return true;
}

bool LGISOlatedLeptonFinderProcessor::IsISOlatedPolynomial( ReconstructedParticle* pfo ) {

	double E     = pfo->getEnergy() ;
	double coneE = getConeEnergy( pfo );

	if ( coneE*coneE <= _isoPolynomialA*E*E + _isoPolynomialB*E + _isoPolynomialC )
		return true ;
	return false;
}

bool LGISOlatedLeptonFinderProcessor::IsISOlatedJet( ReconstructedParticle* pfo ) {
	// jet-based isolated lepton (LAL algorithm)

	if ( _rpJetMap.find( pfo ) == _rpJetMap.end() ) {
		// this is often the case when jet finding fails e.g. due to too few particles in event
		return false;
	}

	ReconstructedParticle* jet = _rpJetMap[pfo];
	TVector3 vec1( pfo->getMomentum() );
	TVector3 jetmom( jet->getMomentum() );
	TLorentzVector jetmom4( jet->getMomentum(), jet->getEnergy() );

	double jetxt = vec1.Pt( jetmom )/jetmom4.M();
	double jetz = pfo->getEnergy()/jet->getEnergy();

	if (jetxt >= _jetIsoVetoMinXt && jetxt < _jetIsoVetoMaxXt
			&& jetz >= _jetIsoVetoMinZ && jetz < _jetIsoVetoMaxZ) {
		//printf("xt=%f z=%f (not pass)\n",jetxt,jetz);
		return false;
	}

	//printf("xt=%f z=%f (PASS)\n",jetxt,jetz);
	return true;
}

bool LGISOlatedLeptonFinderProcessor::PassesImpactParameterCuts( ReconstructedParticle* pfo ) {
	const EVENT::TrackVec & trkvec = pfo->getTracks();

	if (trkvec.size()==0) return false;

	// TODO: more sophisticated pfo/track matching
	double d0 = fabs(trkvec[0]->getD0());
	double z0 = fabs(trkvec[0]->getZ0());
	double r0 = sqrt( d0*d0 + z0*z0 );

	if ( d0 < _minD0 ) return false;
	if ( d0 > _maxD0 ) return false;
	if ( z0 < _minZ0 ) return false;
	if ( z0 > _maxZ0 ) return false;
	if ( r0 < _minR0 ) return false;
	if ( r0 > _maxR0 ) return false;

	return true;
}

bool LGISOlatedLeptonFinderProcessor::PassesImpactParameterSignificanceCuts( ReconstructedParticle* pfo ) {
	const EVENT::TrackVec & trkvec = pfo->getTracks();

	if (trkvec.size()==0) return false;

	// TODO: more sophisticated pfo/track matching
	double d0    = fabs(trkvec[0]->getD0());
	double z0    = fabs(trkvec[0]->getZ0());
	double d0err = sqrt(trkvec[0]->getCovMatrix()[0]);
	double z0err = sqrt(trkvec[0]->getCovMatrix()[9]);

	double d0sig = d0err != 0 ? d0/d0err : 0;
	double z0sig = z0err != 0 ? z0/z0err : 0;
	double r0sig = sqrt( d0sig*d0sig + z0sig*z0sig );

	if ( d0sig < _minD0Sig ) return false;
	if ( d0sig > _maxD0Sig ) return false;
	if ( z0sig < _minZ0Sig ) return false;
	if ( z0sig > _maxZ0Sig ) return false;
	if ( r0sig < _minR0Sig ) return false;
	if ( r0sig > _maxR0Sig ) return false;

	return true;
}

double LGISOlatedLeptonFinderProcessor::getConeEnergy( ReconstructedParticle* pfo ) {
	double coneEC = 0, coneEN=0;

	TVector3 P( pfo->getMomentum() );
	int npfo = _pfoCol->getNumberOfElements();
	for ( int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo_i = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		// don't add itself to the cone energy
		if ( pfo == pfo_i ) continue;

		double charge = pfo_i->getCharge();
		TVector3 P_i( pfo_i->getMomentum() );
		double cosTheta = P.Dot( P_i )/(P.Mag()*P_i.Mag());
		if ( cosTheta >= _cosConeAngle ){
			if( fabs(charge)>0.001)
				coneEC += pfo_i->getEnergy(); 
			else 	
				coneEN += pfo_i->getEnergy();
		}	
	}

	return coneEC+coneEN;
}

void LGISOlatedLeptonFinderProcessor::getCalEnergy( ReconstructedParticle* pfo , double* cale) {
	double ecal = 0;
	double hcal = 0;
	std::vector<lcio::Cluster*> clusters = pfo->getClusters();
	for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
			iCluster!=clusters.end();
			++iCluster) {
		ecal += (*iCluster)->getSubdetectorEnergies()[0];
		hcal += (*iCluster)->getSubdetectorEnergies()[1];
	}
	cale[0] = ecal;
	cale[1] = hcal;
}


