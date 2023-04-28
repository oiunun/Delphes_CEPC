// *******************************************************
// some useful functions
//                ----tianjp
// *******************************************************
#include "Utilities.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <UTIL/LCRelationNavigator.h>

bool  Sort_by_E  (const TLorentzVector& t1, const TLorentzVector& t2) { return t1.E()>t2.E(); }
bool  Sort_by_Ed (const TLorentzVector& t1, const TLorentzVector& t2) { return t1.E()<t2.E(); }

bool  Sort_PFOs_E ( FSParticle* t1,  FSParticle* t2) {
	return (t1->pfo())->getEnergy() > (t2->pfo())->getEnergy(); 
}

bool  Sort_PFOs_Ed ( FSParticle* t1,  FSParticle* t2) {
	return (t1->pfo())->getEnergy() < (t2->pfo())->getEnergy(); 
}

Int_t getMCSerial(MCParticle *mcPart, LCCollection *colMCP) {
	// get the serial code of the particle in the MC particle list
	Int_t nMCP = colMCP->getNumberOfElements();
	Int_t serialID;
	for (Int_t i=0;i<nMCP;i++) {
		MCParticle *mcp = dynamic_cast<MCParticle*>(colMCP->getElementAt(i));
		if (mcp == mcPart) {
			serialID = i;
			return serialID;
		}
	}
	return -1;
}

MCParticle* getLinkedMCParticle(ReconstructedParticle *recPart, LCCollection *colMCTL, Double_t &weight, Int_t &nMCTL) {
	// get the corresponding MC particle of one reconstructed particle using the MCTruthLinker information
	MCParticle *mcLinkedParticle = NULL;
	LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
	LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
	FloatVec vecWgtMCTL = navMCTL->getRelatedToWeights(recPart);
	Double_t mcEnergyMax = -1.0;
	weight = 0.;
	nMCTL = vecMCTL.size();
	for (Int_t i=0;i<nMCTL;i++) {
		MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[i]);
		Double_t mcEnergy = mcPart->getEnergy();
		if (mcEnergy > mcEnergyMax) {   // find the linked particle with largest energy as the mc truth particle
			mcEnergyMax = mcEnergy;
			mcLinkedParticle = mcPart;
			weight = vecWgtMCTL[i];
		}
	}
	return mcLinkedParticle;
}

Int_t getOriginalPDG(MCParticle *mcPart) {
	// get the PDG of the original particle where the MCParticle comes from
	Int_t nParents = mcPart->getParents().size();
	while (nParents > 0) {
		//    MCParticle *mother = mcPart->getParent(0);
		MCParticle *mother = mcPart->getParents()[0];
		//    nParents = mother->getNumberOfParents();
		nParents = mother->getParents().size();
		mcPart = mother;
		Int_t pdg = mcPart->getPDG();
		if (abs(pdg) == 24 || abs(pdg) == 23 || abs(pdg) == 25) break;
	}
	Int_t originalPDG = mcPart->getPDG();
	return originalPDG;
}

Int_t getLeptonID(ReconstructedParticle *recPart) {
	// electron identification using ratios of energies deposited in ECal, HCal and Momentum
	Int_t iLeptonType = 0;
	Double_t fElectronCut1 = 0.5;      // lower edge of totalCalEnergy/momentum
	Double_t fElectronCut2 = 1.3;      // upper edge of totalCalEnergy/momentum
	Double_t fElectronCut3 = 0.9;      // lower edge of ecalEnergy/totalCalEnergy
	Double_t fMuonCut1     = 0.3;      // upper edge of totalCalEnergy/momentum
	//Double_t fMuonCut2   = 0.5;      // upper edge of ecalEnergy/totalCalEnergy
	Double_t fMuonCut3     = -1.2;     // lower edge of yoke energy
	Double_t fPhotonCut1   =  0.7;     // lower edge of totalCalEnergy/momentum
	Double_t fPhotonCut2   =  1.3;     // upper edge of totalCalEnergy/momentum
	Double_t fPhotonCut3   =  0.9;     // lower edge of ecalEnergy/totalCalEnergy

	Double_t charge = recPart->getCharge();
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
	Double_t fEpsilon = 1.E-10;
	if (totalCalEnergy/momentumMagnitude > fElectronCut1 && totalCalEnergy/momentumMagnitude < fElectronCut2 &&
			ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && TMath::Abs(charge) > 0.5) {
		iLeptonType = 11;
	}
	else if (totalCalEnergy/momentumMagnitude > fPhotonCut1 && totalCalEnergy/momentumMagnitude < fPhotonCut2 &&
			ecalEnergy/(totalCalEnergy + fEpsilon) > fPhotonCut3 && TMath::Abs(charge) < 0.5) {
		iLeptonType = 22;
	}
	//  else if (TMath::Abs(charge) > 0.5 && totalCalEnergy/momentumMagnitude < fMuonCut1 && 
	//	   ecalEnergy/(totalCalEnergy + fEpsilon) < fMuonCut2) {
	else if (TMath::Abs(charge) > 0.5 && totalCalEnergy/momentumMagnitude < fMuonCut1 && 
			yokeEnergy > fMuonCut3) {
		iLeptonType = 13;
	}

	return iLeptonType;
}

Bool_t getFSRTag(ReconstructedParticle *motherPart, ReconstructedParticle *recPart, Double_t fCosFSRCut) {
	// tag the particle recPart if it is from the Bremmstrahlung or Final-State-Radiation of another particle motherPart
	Bool_t isFSR = kFALSE;
	//Double_t charge = motherPart->getCharge(); // mother particle should be charged
	TVector3 momentumMother = TVector3(motherPart->getMomentum());
	TVector3 momentum = TVector3(recPart->getMomentum());
	Double_t cosFSR = momentum.Dot(momentumMother)/momentum.Mag()/momentumMother.Mag();
	Int_t iType = getLeptonID(recPart);
	Int_t iTypeMother = getLeptonID(motherPart);
	//  if (TMath::Abs(charge) > 0.5 && cosFSR > fCosFSRCut) {
	//    if (iType == 11 || iType == 22) {
	if (cosFSR > fCosFSRCut && (iTypeMother == 11 || iTypeMother == 13)) {
		if (iType == 22) {
			isFSR = kTRUE;
		}
	}
	return isFSR;  
}

Bool_t getSplitTag(ReconstructedParticle *motherPart, ReconstructedParticle *recPart) {
	// developed from Mark Tomthon's ZFinder Processor
	// tag the particle recPart if it is the cluster splited from motherPart
	Bool_t isSplit = kFALSE;

	// information of mother particle
	Double_t trackMom=0.;
	Double_t clusterMom=0.;
	Double_t sigmaMom=999.;
	Double_t chiMom=999.;
	Double_t sigpMom=0.;
	TVector3 vecMom;
	const EVENT::ClusterVec ci   = motherPart->getClusters();
	const EVENT::TrackVec   ti   = motherPart->getTracks();
	trackMom = motherPart->getEnergy();
	if(ti.size()==1)sigpMom = trackMom*sqrt(ti[0]->getCovMatrix()[5])/TMath::Abs(ti[0]->getOmega());
	if(ci.size()>0){
		clusterMom = ci[0]->getEnergy();
		sigmaMom   = 0.18*sqrt(clusterMom);
		chiMom = (trackMom-clusterMom)/TMath::Sqrt(sigmaMom*sigmaMom+sigpMom*sigpMom);
		vecMom = TVector3(ci[0]->getPosition()[0],ci[0]->getPosition()[1],ci[0]->getPosition()[2]);
	}

	// calculate the distance between cluster and the mother particle cluster
	Double_t dr = 999.;
	const EVENT::ClusterVec c   = recPart->getClusters();
	if(c.size()==1){
		TVector3 vecg(c[0]->getPosition()[0],c[0]->getPosition()[1],c[0]->getPosition()[2]);
		TVector3 v =  vecg.Cross(vecMom);
		float magg = vecg.Mag();
		if(magg>0) dr = v.Mag()/magg;
	}

	// criteria for tag
	Double_t eg = recPart->getEnergy();
	float chiNew = (trackMom-clusterMom-eg)/sigmaMom;
	if(dr<20.0) isSplit = true;
	// if fairly close merge if chi2 for matching improves greatly
	if(dr<30.0 && chiMom>4.0 && TMath::Abs(chiNew)<chiMom) isSplit = true;
	if(dr<40.0 && chiMom>5.0 && TMath::Abs(chiNew)<chiMom) isSplit = true;
	if(dr<50.0 && chiMom>7.0 && TMath::Abs(chiNew)<chiMom) isSplit = true;
	// sanity check
	if(TMath::Abs(chiMom)<2.0 && chiNew*chiNew>chiMom*chiMom+5.0) isSplit = false;
	// always merge if very close - can't expect reconstruction to work
	if(dr<10.0) isSplit = true;

	return isSplit;  
}

Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone) { 
	// get the cone energy of the particle
	Int_t nPFOs = colPFO->getNumberOfElements();
	Double_t coneEnergy = 0.;
	TVector3 momentum0 = TVector3(recPart->getMomentum());
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		if (pfo != recPart) {
			TVector3 momentum = TVector3(pfo->getMomentum());
			Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
			if (cosTheta > cosCone) {
				coneEnergy += pfo->getEnergy();
			}
		}
	}
	return coneEnergy;
}

void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, Double_t coneEnergy[3], Double_t pFSR[4]) { 
	// get the cone energy of the particle
	//  woFSR = kTRUE;
	TLorentzVector lortzFSR = TLorentzVector(0,  0,  0,  0);
	Int_t nPFOs = colPFO->getNumberOfElements();
	TVector3 momentum0 = TVector3(recPart->getMomentum());
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		if (pfo == recPart) continue;
		Double_t energy = pfo->getEnergy();
		TVector3 momentum = TVector3(pfo->getMomentum());
		if (woFSR) {
			Bool_t isFSR = getFSRTag(recPart,pfo);
			if (isFSR) {
				lortzFSR += TLorentzVector(momentum,energy);
				continue;
			}
		}
		Double_t charge = pfo->getCharge();
		Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
		if (cosTheta > cosCone) {
			coneEnergy[0] += energy;
			if (TMath::Abs(charge) < 0.5) {
				coneEnergy[1] += energy;
			}
			else {
				coneEnergy[2] += energy;
			}
		}
	}
	pFSR[0] = lortzFSR.Px();
	pFSR[1] = lortzFSR.Py();
	pFSR[2] = lortzFSR.Pz();
	pFSR[3] = lortzFSR.E();

	return;
}

void getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Bool_t woFSR, Double_t coneEnergy[3], Double_t pFSR[4], 
		Double_t cosCone2, Double_t pCone2[4], Int_t &nConePhoton ) { 
	// get the cone energy of the particle
	// add another larger cone
	//  woFSR = kTRUE;
	TLorentzVector lortzFSR = TLorentzVector(0,  0,  0,  0);
	TLorentzVector lortzCon = TLorentzVector(0,  0,  0,  0);
	Int_t nPFOs = colPFO->getNumberOfElements();
	TVector3 momentum0 = TVector3(recPart->getMomentum());
	nConePhoton = 0;
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		if (pfo == recPart) continue;
		Double_t energy = pfo->getEnergy();
		TVector3 momentum = TVector3(pfo->getMomentum());
		if (woFSR) {
			Bool_t isFSR = getFSRTag(recPart,pfo);
			if (isFSR) {
				lortzFSR += TLorentzVector(momentum,energy);
				Int_t iType = getLeptonID(pfo);
				if (iType == 22) nConePhoton++;
				continue;
			}
		}
		Double_t charge = pfo->getCharge();
		Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
		if (cosTheta > cosCone) {
			coneEnergy[0] += energy;
			if (TMath::Abs(charge) < 0.5) {
				coneEnergy[1] += energy;
			}
			else {
				coneEnergy[2] += energy;
			}
			Int_t iType = getLeptonID(pfo);
			if (iType == 22) nConePhoton++;
		}
		if (cosTheta > cosCone2) {
			lortzCon += TLorentzVector(momentum,energy);
		}
	}
	pFSR[0] = lortzFSR.Px();
	pFSR[1] = lortzFSR.Py();
	pFSR[2] = lortzFSR.Pz();
	pFSR[3] = lortzFSR.E();
	pCone2[0] = lortzCon.Px();
	pCone2[1] = lortzCon.Py();
	pCone2[2] = lortzCon.Pz();
	pCone2[3] = lortzCon.E();
	return;
}

TLorentzVector getFSRMomentum(ReconstructedParticle *recPart, LCCollection *colPFO) { 
	// get the cone energy of the particle
	TLorentzVector lortzFSR = TLorentzVector(0.,0.,0.,0.);
	Int_t nPFOs = colPFO->getNumberOfElements();
	TVector3 momentum0 = TVector3(recPart->getMomentum());
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		if (pfo == recPart) continue;
		Bool_t isFSR = getFSRTag(recPart,pfo);
		if (! isFSR) continue;
		Double_t energy = pfo->getEnergy();
		TVector3 momentum = TVector3(pfo->getMomentum());
		lortzFSR += TLorentzVector(momentum,energy);
	}

	return lortzFSR;
}

Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, Int_t mode) { 
	// get the cone energy of the particle
	Int_t nPFOs = colPFO->getNumberOfElements();
	Double_t coneEnergy = 0.,coneEnergyC = 0.,coneEnergyN = 0.;
	TVector3 momentum0 = TVector3(recPart->getMomentum());
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		if (pfo != recPart) {
			TVector3 momentum = TVector3(pfo->getMomentum());
			Int_t iCharge = (Int_t)pfo->getCharge();
			Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
			if (cosTheta > cosCone) {
				coneEnergy += pfo->getEnergy();
				if (iCharge == 0) {
					coneEnergyN += pfo->getEnergy();
				}
				else {
					coneEnergyC += pfo->getEnergy();
				}
			}
		}
	}
	if (mode == 0) {
		return coneEnergy;
	}
	else if (mode == 1) {
		return coneEnergyC;
	}
	else if (mode == 2) {
		return coneEnergyN;
	}
	else {
		return 99999.;
	}
}

Double_t getConeEnergy(ReconstructedParticle *recPart, LCCollection *colPFO, Double_t cosCone, 
		std::vector<lcio::ReconstructedParticle*> &conePFOs) {
	// get the cone energy of the particle
	Int_t nPFOs = colPFO->getNumberOfElements();
	Double_t coneEnergy = 0.;
	TVector3 momentum0 = TVector3(recPart->getMomentum());
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		if (pfo != recPart) {
			TVector3 momentum = TVector3(pfo->getMomentum());
			Double_t cosTheta = momentum0.Dot(momentum)/momentum0.Mag()/momentum.Mag();
			if (cosTheta > cosCone) {
				coneEnergy += pfo->getEnergy();
				conePFOs.push_back(pfo);
			}
		}
	}
	return coneEnergy;
}

Double_t getInvariantMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2) {
	// get the invariant mass of two particles
	TVector3 momentum1 = TVector3(recPart1->getMomentum());
	TVector3 momentum2 = TVector3(recPart2->getMomentum());
	Double_t energy1 = recPart1->getEnergy();
	Double_t energy2 = recPart2->getEnergy();
	Double_t invariantMass = 0.;
	Double_t invariantMass2 = (energy1+energy2)*(energy1+energy2)-(momentum1+momentum2).Mag2();
	if (invariantMass2 > 0.) {
		invariantMass = sqrt(invariantMass2);
	}
	else {
		invariantMass = -sqrt(TMath::Abs(invariantMass2));
	}
	return invariantMass;
}

Double_t        getRecoilingMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2, Double_t ecms){

	// get the recoil mass of two particles
	TVector3 momentum1 = TVector3(recPart1->getMomentum());
	TVector3 momentum2 = TVector3(recPart2->getMomentum());
	Double_t energy1 = recPart1->getEnergy();
	Double_t energy2 = recPart2->getEnergy();
	Double_t invariantMass = 0.;
	Double_t invariantMass2 = (ecms-energy1-energy2)*(ecms-energy1-energy2)-(momentum1+momentum2).Mag2();
	if (invariantMass2 > 0.) {
		invariantMass = sqrt(invariantMass2);
	}
	else {
		invariantMass = -sqrt(TMath::Abs(invariantMass2));
	}
	return invariantMass;

}
Double_t getInvariantMass(ReconstructedParticle *recPart1, ReconstructedParticle *recPart2, ReconstructedParticle *recPart3) {
	// get the invariant mass of three particles
	TVector3 momentum1 = TVector3(recPart1->getMomentum());
	TVector3 momentum2 = TVector3(recPart2->getMomentum());
	TVector3 momentum3 = TVector3(recPart3->getMomentum());
	Double_t energy1 = recPart1->getEnergy();
	Double_t energy2 = recPart2->getEnergy();
	Double_t energy3 = recPart3->getEnergy();
	Double_t invariantMass = 0.;
	Double_t invariantMass2 = (energy1+energy2+energy3)*(energy1+energy2+energy3)-(momentum1+momentum2+momentum3).Mag2();
	if (invariantMass2 > 0.) {
		invariantMass = sqrt(invariantMass2);
	}
	else {
		invariantMass = -sqrt(TMath::Abs(invariantMass2));
	}
	return invariantMass;
}

Int_t isSelectedByFastJet( ReconstructedParticle *pfo, LCCollection *colFastJet, Double_t &ratioEPartEJet, Double_t &ratioPTMJet ) { 
	// check the PFO if it is selectred by FastJet clustering
	Int_t iFastJet = 0;
	Int_t nJets = colFastJet->getNumberOfElements();
	Int_t iJet = -1;
	for (Int_t i=0;i<nJets;i++) {
		ReconstructedParticle *jet = dynamic_cast<ReconstructedParticle*>(colFastJet->getElementAt(i));
		std::vector<lcio::ReconstructedParticle*> partVec = jet->getParticles();
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPart=partVec.begin();iPart!=partVec.end();++iPart) {
			if ((*iPart) == pfo) {
				iFastJet = 1;
				iJet = i;
				break;
			}
		}
	}
	if (iJet >= 0) {
		ReconstructedParticle *theJet = dynamic_cast<ReconstructedParticle*>(colFastJet->getElementAt(iJet)); // the jet where the pfo belongs
		// get the variables used by LAL Lepton Finder
		Double_t ePart = pfo->getEnergy();
		TVector3 pPart = pfo->getMomentum();
		Double_t eJet  = theJet->getEnergy();
		TVector3 pJet  = theJet->getMomentum();
		//    Double_t mJet  = theJet->getMass();
		Double_t mJet  = eJet*eJet > pJet.Mag2() ? TMath::Sqrt(eJet*eJet-pJet.Mag2()) : TMath::Sqrt(-eJet*eJet+pJet.Mag2());
		ratioEPartEJet = ePart/eJet;
		ratioPTMJet = pPart.Pt(pJet)/mJet;
	}

	return iFastJet;
}

Double_t getEnergyWeightedIP( LCCollection *colPFO ) {
	// get the average ip position weighted by energy of each particle

	Double_t z0IP = 0.;
	Double_t sumWz0 = 0., sumE = 0.;
	Int_t nPFOs = colPFO->getNumberOfElements();
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		Double_t energy = recPart->getEnergy();
		//Double_t charge = recPart->getCharge();
		TrackVec tckvec = recPart->getTracks();
		Int_t ntracks = tckvec.size();
		if (ntracks > 0) {
			Double_t z0 = tckvec[0]->getZ0();
			sumWz0 += z0*energy;
			sumE += energy;
		}
	}  
	if (sumE > 0.) {
		z0IP = sumWz0/sumE;
	}
	else {
		z0IP = 0.;
	}

	return z0IP;
}
