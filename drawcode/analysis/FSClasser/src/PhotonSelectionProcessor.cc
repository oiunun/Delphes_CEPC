// *****************************************************
// e+e- ------> ZH ------> (qq)(gamma gamma)
// Processor for gamma selection
//                        ----LIGang
// *****************************************************
#include "PhotonSelectionProcessor.h"
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

#include "Utilities.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

PhotonSelectionProcessor aPhotonSelectionProcessor ;


PhotonSelectionProcessor::PhotonSelectionProcessor() : Processor("PhotonSelectionProcessor") {

	// modify processor description
	_description = "PhotonSelectionProcessor does whatever it does ..." ;


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

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"OutputNewPFOsCollection" , 
			"Name of the NewPFOs collection"  ,
			_colNewPFOs ,
			std::string("NewPFOs") ) ;

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputLeptonsCollection",
			"Name of collection with the selected leptons",
			_colLeptons,
			std::string("leptons") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputPhotonsCollection",
			"Name of collection with the selected photons",
			_colPhotons,
			std::string("photons") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputZPFOsCollection",
			"Name of collection with the selected Z PFOs",
			_colZPFOs,
			std::string("zpfos") );
}

void PhotonSelectionProcessor::init() { 

	streamlog_out(DEBUG) << "   init called  " 
		<< std::endl ;


	printParameters() ;

	_nRun = 0 ;
	_nEvt = 0 ;

	hStat = 0;
	//  TFile *outRootFile = new TFile("output.root","RECREATE");

}

void PhotonSelectionProcessor::processRunHeader( LCRunHeader* run) { 

	_nRun++ ;
} 

void PhotonSelectionProcessor::processEvent( LCEvent * evt ) { 


	// this gets called for every event 
	// usually the working horse ...
	_nEvt++;

#if 1
	Double_t fEtrackCut = -1.;            // lower edge of each PFO energy 
#else
	Double_t fEtrackCut = 0.05;           // lower edge of each PFO energy 
#endif
	Double_t fElectronCut1 = 0.5;         // lower edge of totalCalEnergy/momentum
	Double_t fElectronCut2 = 1.3;         // upper edge of totalCalEnergy/momentum
	Double_t fElectronCut3 = 0.8; //0.9;         // lower edge of ecalEnergy/totalCalEnergy
	Double_t fMuonCut1 = 0.3;             // upper edge of totalCalEnergy/momentum
	Double_t fMuonCut3 = 1.2;             // lower edge of yoke energy
	Double_t fCosConeCut = 0.98;          // the angle of cone around the direction of pfo
	Double_t kMassZ = 91.187;             // Z mass
	Double_t kMassH = 125.;               // Higgs mass
	Double_t kMassZOff = kMassH - kMassZ; // Z mass
	Double_t fMassZCut = 40.;             // mass cut for lepton pair from Z
	Double_t fMassHCut = 30.;             // mass cut for photon pair from H
	Double_t fEpsilon = 1.E-10;
	// Fisher coefficients
	Double_t c0Electron = 12.2;           //used in llHH mode
	Double_t c1Electron = 0.87; 
	Double_t c0Muon = 12.6; 
	Double_t c1Muon = 4.62; 

	//cerr << endl << "Hello, Photon Selection!" << endl;

	// -- Get the MCTruth Linker --
	LCCollection *colMCTL = evt->getCollection(_colMCTL);
	LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

	// -- Read out PFO information --
	LCCollection *colPFO = evt->getCollection(_colPFOs);
	if (!colPFO) {
		std::cerr << "No PFO Collection Found!" << std::endl;
		throw marlin::SkipEventException(this);
	}
	Int_t nPFOs = colPFO->getNumberOfElements();
	//cerr << "Number of PFOs: " << nPFOs << endl;
	LCCollectionVec *pNewPFOsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollectionVec *pLeptonsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollectionVec *pPhotonsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	pNewPFOsCollection->setSubset(true);
	pLeptonsCollection->setSubset(true);
	pPhotonsCollection->setSubset(true);
	std::vector<lcio::ReconstructedParticle*> newPFOs;
	std::vector<lcio::ReconstructedParticle*> electrons;
	std::vector<lcio::ReconstructedParticle*> muons;
	std::vector<lcio::ReconstructedParticle*> photons;


	//  Bool_t iLepTune = kTRUE;
	Bool_t iLepTune = kFALSE;
	// loop all the PFOs
	Double_t energyLep1=0.,energyLep2=0.;
	Double_t energyLep1Fsr=0.,energyLep2Fsr=0.;
	Double_t cosLep1=0.,cosLep2=0.;
	Double_t nLep1Found=0.,nLep2Found=0.;
	TVector3 momentumLep1 = TVector3(0.,0.,0.);
	TVector3 momentumLep2 = TVector3(0.,0.,0.);
	Int_t indexLep1=-1,indexLep2=-1;
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
		FloatVec vecWgtMCTL = navMCTL->getRelatedToWeights(recPart);
		Int_t mcpdg,motherpdg,mmotherpdg;
		Double_t mcwgt=0.;
		mcpdg = 0;
		motherpdg = -99999;
		mmotherpdg = -99999;
		Double_t deltaE = -99999.;
		Double_t energyLink = -99999.;
		Int_t mcoriginal = 0;
		Int_t mcndaughters = 0;
		Int_t nMCTL = vecMCTL.size();
		Int_t iOverlay = 0;
		if (vecMCTL.size() > 0) {
			MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
			if (mcPart->isOverlay()) iOverlay = 1;
			mcpdg = mcPart->getPDG();
			mcwgt = vecWgtMCTL[0];
			deltaE = mcPart->getEnergy()-recPart->getEnergy();
			energyLink = mcPart->getEnergy();
			mcoriginal = getOriginalPDG(mcPart);
			motherpdg  = 0;
			mcndaughters = mcPart->getDaughters().size();
			if (mcPart->getParents().size() != 0) {
				MCParticle *motherPart = mcPart->getParents()[0];
				motherpdg = motherPart->getPDG();
				mmotherpdg = 0;
				if (motherPart->getParents().size() != 0) {
					MCParticle *mmotherPart = motherPart->getParents()[0];
					mmotherpdg = mmotherPart->getPDG();
				}
			}
		}
		Double_t energy = recPart->getEnergy();
		Double_t charge = recPart->getCharge();
		Int_t itype = recPart->getType();
		Int_t pid = 0;
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
		Double_t r0 = TMath::Sqrt(d0*d0+z0*z0);
		Double_t nsigr0 = TMath::Sqrt(nsigd0*nsigd0+nsigz0*nsigz0);
		
		if (energy > fEtrackCut) {
			newPFOs.push_back(recPart);
			Double_t ecalEnergy = 0;
			Double_t hcalEnergy = 0;
			Double_t yokeEnergy = 0;
			Double_t totalCalEnergy = 0;
			Int_t nHits = 0;
			std::vector<lcio::Cluster*> clusters = recPart->getClusters();
			for (std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster) {
				ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
				hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
				yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
				ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
				hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
				CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
				nHits += calHits.size();
			}
			totalCalEnergy = ecalEnergy + hcalEnergy;
			TVector3 momentum = TVector3(recPart->getMomentum());
			Double_t momentumMagnitude = momentum.Mag();
			Double_t cosTheta = momentum.CosTheta();
			//get cone information
			Bool_t woFSR = kTRUE;
			Double_t coneEnergy0[3] = {0.,0.,0.};
			Double_t pFSR[4] = {0.,0.,0.,0.};
			getConeEnergy(recPart,colPFO,fCosConeCut,woFSR,coneEnergy0,pFSR);
			Double_t coneEnergy = coneEnergy0[0];
			Double_t coneEN     = coneEnergy0[1];
			Double_t coneEC     = coneEnergy0[2];
			TLorentzVector lortzFSR = TLorentzVector(pFSR[0],pFSR[1],pFSR[2],pFSR[3]);
			Double_t energyCorr = energy + lortzFSR.E();
			TVector3 momentumCorr = momentum + TVector3(lortzFSR.Px(),lortzFSR.Py(),lortzFSR.Pz());
			Double_t momentumMagCorr = momentumCorr.Mag();

			Int_t nConePFOs    = 0;
			Int_t nConeCharged = 0;
			Int_t nConeNeutral = 0;

			if (charge == -1 && abs(mcoriginal) == 11 && (mcpdg == +11 || mcpdg == +13) && motherpdg!=22 && energy > energyLep1) {
				energyLep1 = energy;
				energyLep1Fsr = energyCorr;
				cosLep1    = cosTheta;
				nLep1Found++;
				indexLep1 = i;
				momentumLep1 = momentum;
			}
			if (charge == +1 && abs(mcoriginal) == 11 && (mcpdg == -11 || mcpdg == -13) && motherpdg!=22 && energy > energyLep2) {
				energyLep2 = energy;
				energyLep2Fsr = energyCorr;
				cosLep2    = cosTheta;
				nLep2Found++;
				indexLep2 = i;
				momentumLep2 = momentum;
			}
			// select the leptons
			if (charge != 0 && 
					totalCalEnergy/momentumMagnitude > fElectronCut1 && totalCalEnergy/momentumMagnitude < fElectronCut2 &&
					ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && 
					(momentumMagnitude > c0Electron + c1Electron*coneEC)) {
				if (nsigd0 > 50 || nsigz0 > 5) continue;   // contraint to primary vertex
				electrons.push_back(recPart);
			}
			if (charge != 0 && 
					totalCalEnergy/momentumMagnitude < fMuonCut1 && 
					yokeEnergy > fMuonCut3 && 
					(momentumMagnitude > c0Muon + c1Muon*coneEC)) {
				if (nsigd0 > 5 || nsigz0 > 5) continue;  // contraint to primary vertex
				muons.push_back(recPart);
			}
			// select the photons
			if (charge == 0 && 
					ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && totalCalEnergy > fElectronCut2 
				) {
				photons.push_back(recPart);
			}
		}
	}
	Int_t nelectrons = electrons.size();
	Int_t nmuons     = muons.size();
	Int_t nphotons   = photons.size();
	//  cerr << "nelectrons: " << nelectrons << "  nmuons: " << nmuons << endl;
	Int_t iLeptonType = 0;
	if (nelectrons > 0) iLeptonType += 1;
	if (nmuons > 0) iLeptonType += 10;
	//  if (nelectrons < 2 && nmuons < 2) throw marlin::SkipEventException(this);
	//  if (nelectrons+nmuons > 2) throw marlin::SkipEventException(this);

	// find the lepton pair nearest to the Z mass
	std::vector<lcio::ReconstructedParticle*> electronZ;
	std::vector<lcio::ReconstructedParticle*> muonZ;
	std::vector<lcio::ReconstructedParticle*> leptonZ;
	std::vector<lcio::ReconstructedParticle*> photonH;
	std::vector<lcio::ReconstructedParticle*> photonSplit;
	Double_t deltaMassZ = fMassZCut;
	Double_t massZ=0.;
	if (electrons.size() > 1) {
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iElectron=electrons.begin();iElectron<electrons.end()-1;iElectron++) {
			for (std::vector<lcio::ReconstructedParticle*>::const_iterator jElectron=iElectron+1;jElectron<electrons.end();jElectron++) {
				if ((*iElectron)->getCharge() != (*jElectron)->getCharge()) {
					Double_t mass = getInvariantMass((*iElectron),(*jElectron));
					if (TMath::Abs(mass-kMassZ) < deltaMassZ+20 || TMath::Abs(mass-kMassZOff) < deltaMassZ+20) {
						deltaMassZ = TMath::Abs(mass-kMassZ) < TMath::Abs(mass-kMassZOff) ? TMath::Abs(mass-kMassZ) : TMath::Abs(mass-kMassZOff);
						massZ = mass;
						electronZ.clear();
						electronZ.push_back(*iElectron);
						electronZ.push_back(*jElectron);
					}
				}
			}
		}
	}
	if (muons.size() > 1) {
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iMuon=muons.begin();iMuon<muons.end()-1;iMuon++) {
			for (std::vector<lcio::ReconstructedParticle*>::const_iterator jMuon=iMuon+1;jMuon<muons.end();jMuon++) {
				if ((*iMuon)->getCharge() != (*jMuon)->getCharge()) {
					Double_t mass = getInvariantMass((*iMuon),(*jMuon));
					if (TMath::Abs(mass-kMassZ) < deltaMassZ || TMath::Abs(mass-kMassZOff) < deltaMassZ) {
						deltaMassZ = TMath::Abs(mass-kMassZ) < TMath::Abs(mass-kMassZOff) ? TMath::Abs(mass-kMassZ) : TMath::Abs(mass-kMassZOff);
						massZ = mass;
						muonZ.clear();
						muonZ.push_back(*iMuon);
						muonZ.push_back(*jMuon);
					}
				}
			}
		}
	}

	Double_t deltaMassH = fMassHCut;
	Double_t massH=0.;
	if (photons.size() > 1) {
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPhoton=photons.begin();iPhoton<photons.end()-1;iPhoton++) {
			for (std::vector<lcio::ReconstructedParticle*>::const_iterator jPhoton=iPhoton+1;jPhoton<photons.end();jPhoton++) {
				Double_t mass = getInvariantMass((*iPhoton),(*jPhoton));
				if ( TMath::Abs(mass-kMassH) < deltaMassH ) {
					deltaMassH = TMath::Abs(mass-kMassH);
					massH = mass;
					photonH.clear();
					photonH.push_back(*iPhoton);
					photonH.push_back(*jPhoton);
				}
			}
		}
	}
	//
	//cout<<"No of photons = "<<photonH.size()<<endl;
	//
	Double_t energyLep1Rec=0.,energyLep2Rec=0.;
	Double_t energyLep1FSR=0.,energyLep2FSR=0.;
	Int_t mcOriginalLep1=0,mcOriginalLep2=0;
	if (muonZ.size() == 0 && electronZ.size() == 0 && photonH.size()==0) {
		//cerr << "no lepton pair candidate found!" << endl;
		//throw marlin::SkipEventException(this);
	} else if (muonZ.size() > 0) {
		iLeptonType = 2;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iMuon=muonZ.begin();iMuon<muonZ.end();iMuon++) {
			leptonZ.push_back(*iMuon);
			Int_t mcoriginal = 0;
			LCObjectVec vecMCTL = navMCTL->getRelatedToObjects((*iMuon));
			//Int_t nMCTL = vecMCTL.size();
			if (vecMCTL.size() > 0) {
				MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
				mcoriginal = getOriginalPDG(mcPart);
			}
			if ((*iMuon)->getCharge() < 0) {
				energyLep1FSR += (*iMuon)->getEnergy();
				energyLep1Rec += (*iMuon)->getEnergy();
				mcOriginalLep1 = mcoriginal;
			}
			else {
				energyLep2FSR += (*iMuon)->getEnergy();
				energyLep2Rec += (*iMuon)->getEnergy();
				mcOriginalLep2 = mcoriginal;
			}
			for (Int_t i=0;i<nPFOs;i++) { // recover the FSR for muon
				ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
				if (recPart == (*iMuon)) continue;
				Bool_t isLep = kFALSE;
				for (std::vector<lcio::ReconstructedParticle*>::const_iterator iLep=leptonZ.begin();iLep<leptonZ.end();++iLep) {
					if (recPart == (*iLep)) isLep = kTRUE;
				}
				if (isLep) continue;
				Bool_t isFSR = getFSRTag((*iMuon),recPart);
				if (! isFSR) continue;
				leptonZ.push_back(recPart);
				if ((*iMuon)->getCharge() < 0) {
					energyLep1FSR += recPart->getEnergy();
				}
				else {
					energyLep2FSR += recPart->getEnergy();
				}
			}
			pLeptonsCollection->addElement(*iMuon);
		}
	} else if (electronZ.size() > 0) {
		iLeptonType = 1;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iElectron=electronZ.begin();iElectron<electronZ.end();iElectron++) {
			leptonZ.push_back(*iElectron);
			Int_t mcoriginal = 0;
			LCObjectVec vecMCTL = navMCTL->getRelatedToObjects((*iElectron));
			//Int_t nMCTL = vecMCTL.size();
			if (vecMCTL.size() > 0) {
				MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[0]);
				mcoriginal = getOriginalPDG(mcPart);
			}
			if ((*iElectron)->getCharge() < 0) {
				energyLep1FSR += (*iElectron)->getEnergy();
				energyLep1Rec += (*iElectron)->getEnergy();
				mcOriginalLep1 = mcoriginal;
			}
			else {
				energyLep2FSR += (*iElectron)->getEnergy();
				energyLep2Rec += (*iElectron)->getEnergy();
				mcOriginalLep2 = mcoriginal;
			}
			for (Int_t i=0;i<nPFOs;i++) {  // recover the FSR for electron
				ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
				if (recPart == (*iElectron)) continue;
				Bool_t isLep = kFALSE;
				for (std::vector<lcio::ReconstructedParticle*>::const_iterator iLep=leptonZ.begin();iLep<leptonZ.end();++iLep) {
					if (recPart == (*iLep)) isLep = kTRUE;
				}
				if (isLep) continue;
				Bool_t isFSR = getFSRTag((*iElectron),recPart);
				if (! isFSR) continue;
				leptonZ.push_back(recPart);
				Bool_t isSplit = getSplitTag((*iElectron),recPart);
				if (isSplit) photonSplit.push_back(recPart);
				if ((*iElectron)->getCharge() < 0) {
					energyLep1FSR += recPart->getEnergy();
				}
				else {
					energyLep2FSR += recPart->getEnergy();
				}
			}
			pLeptonsCollection->addElement(*iElectron);
		}
	}
	else if ( photonH.size()>0){
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPhoton=photonH.begin();iPhoton<photonH.end();iPhoton++) {
			pPhotonsCollection->addElement(*iPhoton);
		}
	}

	// create the ZPFOs collection
	LCCollectionVec * pZCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	ReconstructedParticleImpl * recoZ = new ReconstructedParticleImpl();
	TLorentzVector lortzZ_FSR = TLorentzVector(0.,0.,0.,0.);
	TLorentzVector lortzZ_FSR_Split = TLorentzVector(0.,0.,0.,0.);
	for (std::vector<lcio::ReconstructedParticle*>::const_iterator iLep=leptonZ.begin();iLep<leptonZ.end();++iLep) {
		lortzZ_FSR_Split += TLorentzVector((*iLep)->getMomentum(),(*iLep)->getEnergy());
		// do not add the energy and momentum if it is splitted cluster
		Bool_t isSplit = kFALSE;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iSp=photonSplit.begin();iSp<photonSplit.end();++iSp) {
			if ((*iLep) == (*iSp)) isSplit = kTRUE;
		}
		if (isSplit) continue;
		recoZ->addParticle(*iLep);
		lortzZ_FSR += TLorentzVector((*iLep)->getMomentum(),(*iLep)->getEnergy());
	}
	Double_t energyZ_FSR = lortzZ_FSR.E();
	Double_t massZ_FSR = lortzZ_FSR.M();
	Double_t momentumZ_FSR[3] = {lortzZ_FSR.Px(),lortzZ_FSR.Py(),lortzZ_FSR.Pz()};
	recoZ->setMomentum(momentumZ_FSR);
	recoZ->setEnergy(energyZ_FSR);
	recoZ->setMass(massZ_FSR);
	recoZ->setCharge(0.);
	recoZ->setType(94);
	pZCollection->addElement(recoZ);
	Double_t massZ_Split = lortzZ_FSR_Split.M();
	Int_t nSplit = photonSplit.size();

	// save the other PFOs to a new collection
	for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=newPFOs.begin();iObj<newPFOs.end();++iObj) {
		Bool_t isPho=kFALSE;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPho=photonH.begin();iPho<photonH.end();++iPho) {
			if ((*iObj) == (*iPho)) isPho = kTRUE;
		}
		if (!isPho) pNewPFOsCollection->addElement(*iObj);
	}
	evt->addCollection(pLeptonsCollection, _colLeptons.c_str());
	evt->addCollection(pZCollection,       _colZPFOs.c_str()  );

	evt->addCollection(pNewPFOsCollection, _colNewPFOs.c_str());
	evt->addCollection(pPhotonsCollection, _colPhotons.c_str());

	Int_t nNewPFOs = pNewPFOsCollection->getNumberOfElements();

	//-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

	_nEvt ++ ;

}



void PhotonSelectionProcessor::check( LCEvent * evt ) { 
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void PhotonSelectionProcessor::end(){ 

	cerr << "PhotonSelectionProcessor::end()  " << name() 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< endl ;
	//  cerr << endl;
	cerr << "  =============" << endl;
	cerr << "   Cut Summary " << endl;
	cerr << "  =============" << endl;
	cerr << "   ll+4 Jet    " << endl;
	cerr << "  =============" << endl;
	cerr << endl;
}
