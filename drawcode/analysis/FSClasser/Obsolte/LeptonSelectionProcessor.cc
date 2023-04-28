// *****************************************************
// e+e- ------> vvH ------> (vv)(WW*)------> (vv)(qqqq)
// Processor for leptons selection
//                        ----Junping
// *****************************************************
#include "LeptonSelectionProcessor.h"
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

LeptonSelectionProcessor aLeptonSelectionProcessor ;


LeptonSelectionProcessor::LeptonSelectionProcessor() : Processor("LeptonSelectionProcessor") {

	// modify processor description
	_description = "LeptonSelectionProcessor does whatever it does ..." ;


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
			"OutputIsoLeptonsCollection",
			"Name of collection with the selected isolated leptons",
			_colIsoLeptons,
			std::string("IsoLeptons") );


	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputZPFOsCollection",
			"Name of collection with the selected Z PFOs",
			_colZPFOs,
			std::string("zpfos") );
}

void LeptonSelectionProcessor::init() { 

	streamlog_out(DEBUG) << "   init called  " 
		<< std::endl ;


	// usually a good idea to
	printParameters() ;

	_nRun = 0 ;
	_nEvt = 0 ;

	hStat = 0;
	//  TFile *outRootFile = new TFile("output.root","RECREATE");

}

void LeptonSelectionProcessor::processRunHeader( LCRunHeader* run) { 

	_nRun++ ;
} 

void LeptonSelectionProcessor::processEvent( LCEvent * evt ) { 


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
	Double_t fElectronCut3 = 0.9;         // lower edge of ecalEnergy/totalCalEnergy
	Double_t fMuonCut1     = 0.3;         // upper edge of totalCalEnergy/momentum
	Double_t fMuonCut3     =-1.2;         // lower edge of yoke energy
	Double_t fCosConeCut   =0.98;         // the angle of cone around the direction of pfo
	// Fisher coefficients
	Double_t c0Electron    = 12.2;        //used in llHH mode
	Double_t c1Electron    = 0.87; 
	Double_t c0Muon        = 12.6; 
	Double_t c1Muon        = 4.62;
	//
	Double_t kMassZ    = 91.187;          // Z mass
	Double_t kMassH    = 125.;            // Higgs mass
	Double_t kMassZOff = kMassH - kMassZ; // Z mass
	Double_t fMassZCut = 40.;             // mass cut for lepton pair from Z
	Double_t fEpsilon  = 1.E-10;
	//
	//	
	// cut table
	if (!hStat) hStat = new TH1D("hStat", "Cut Table", 20, 0, 20);
	Double_t selid = -0.5;
	hStat->Fill(++selid);
	gCutName[(Int_t)selid] << "No Cuts" << ends;

	TDirectory *last = gDirectory;
	gFile->cd("/");

	//cerr << endl << "Hello, Lepton Selection!" << endl;

	static TNtupleD *hEvt = 0;
	if (!hEvt) {
		cerr << "First Event!" << endl;
		stringstream tupstr;
		tupstr 
			<< "nMCP:nPFO:nelectron:nmuon:ltype"                << ":"
			<< "nhbb:evis:ptot:px:py:pz:pt"                     << ":"
			<< "nhww:nzqq:nzee:nzmm:nn1:nn2:nn3"                << ":"
			<< "nzvv:nhgg:nhcc:nhzz:nhtt"                       << ":"
			<< "zmass:zmassr"                                   << ":"
			<< "elep1:elep2:elepf1:elepf2"                      << ":"
			<< "elepr1:elepr2"                                  << ":"
			<< "elep1mc:elep2mc:cos1mc:cos2mc:mcorig1:mcorig2"  << ":"
			<< "zmassf:nsplit:npfonew"
			<< ends;
		hEvt = new TNtupleD("hEvt","",tupstr.str().data());
	}

	// -- Get the MCTruth Linker --
	LCCollection *colMCTL = evt->getCollection(_colMCTL);
	LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

	// -- Read out MC information --  
	LCCollection *colMC = evt->getCollection(_colMCP);
	if (!colMC) {
		std::cerr << "No MC Collection Found!" << std::endl;
		throw marlin::SkipEventException(this);
	}
	Int_t nMCP = colMC->getNumberOfElements();
	static TNtupleD *hMc = 0;
	if (!hMc) {
		stringstream tupstr_mc;
		tupstr_mc << "id:pdg:motherpdg:mass:nparents:energy:ndaughters:original:daughterpdg" << ":"
			<< "charge"
			<< ends;
		hMc = new TNtupleD("hMc","",tupstr_mc.str().data());
	}
	static TNtupleD *hGen = 0;
	if (!hGen) {
		stringstream tupstr_gen;
		tupstr_gen << "nn1:nn2:nn3"         << ":"
			<< "nhbb:nhww:nzqq:nzvv:nzee:nzmm"   << ":"
			<< "nhgg:nhcc:nhzz:nhtt"
			<< ends;
		hGen = new TNtupleD("hGen","",tupstr_gen.str().data());
	}
	Double_t elepton1=0.,elepton2=0.,elepton3=0.,elepton4=0.;
	Double_t energyLep1MC=0.,energyLep2MC=0.;
	Double_t cosLep1MC=0.,cosLep2MC=0.;
	Double_t cos1=0.,cos2=0.,cos3=0.,cos4=0.;
	Double_t phi1=0.,phi2=0.,phi3=0.,phi4=0.;
	Int_t    nMCTL1=0,nMCTL2=0,nMCTL3=0,nMCTL4=0;
	Double_t chargeRec1=99.,chargeRec2=99.,chargeRec3=99.,chargeRec4=99.;
	Double_t wgtRec1=0.,wgtRec2=0.,wgtRec3=0.,wgtRec4=0.;
	Int_t    nHbb = 0;  // tag H---> b b
	Int_t    nHww = 0;  // tag H---> W W*
	Int_t    nZqq = 0, nZvv = 0, nZee = 0, nZmm = 0; // tag Z decay
	Int_t    nHtt = 0;  // tag H---> tau tau
	Int_t    nHgg = 0;
	Int_t    nHzz = 0;
	Int_t    nHcc = 0;
	Int_t    nn1=0,nn2=0,nn3=0; // tag fusion or strahlung
	TLorentzVector lortzLep1,lortzLep2,lortzZ;
	Bool_t iMCDebug = kFALSE;
	//  Bool_t iMCDebug = kTRUE;
	if (iMCDebug) {
		cerr << "Serial" << "    " << "PDG" << " " << "Mother" << "    " << "Charge" << "      " << "Mass" << "     " << "Energy" << " " 
			<< "NumberOfDaughters" << " " << "NumberOfParents" << "  " << "Original" << endl;
	}
	for (Int_t i=0;i<nMCP;i++) {
		MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
		Int_t pdg = mcPart->getPDG();
		Int_t nparents = mcPart->getParents().size();
		Int_t motherpdg = 0, mmotherpdg = 0;
		if (nparents > 0) {
			MCParticle *mother = mcPart->getParents()[0];
			motherpdg = mother->getPDG();
			if (mother->getParents().size() > 0) {
				MCParticle *mmother = mother->getParents()[0];
				mmotherpdg = mmother->getPDG();
			}
		}
		Double_t charge = mcPart->getCharge();
		Double_t mass = mcPart->getMass();
		Double_t energy = mcPart->getEnergy();
		TVector3 pv = TVector3(mcPart->getMomentum());
		Int_t ndaughters = mcPart->getDaughters().size();
		Int_t daughterpdg = 0;
		if (ndaughters > 0) {
			MCParticle *daughter = mcPart->getDaughters()[0];
			daughterpdg = daughter->getPDG();
		}
		TLorentzVector lortz = TLorentzVector(pv,energy);
		Int_t originalPDG = getOriginalPDG(mcPart);
		Double_t datamc[100];
		datamc[0] = i;
		datamc[1] = pdg;
		datamc[2] = motherpdg;
		datamc[3] = mass;
		datamc[4] = nparents;
		datamc[5] = energy;
		datamc[6] = ndaughters;
		datamc[7] = originalPDG;
		datamc[8] = daughterpdg;
		datamc[9] = charge;
		//    hMc->Fill(datamc);
		if (iMCDebug) {
			cerr << setw(6) << i << setw(7) << pdg << setw(7) << motherpdg << setw(10) << charge << setw(10) << mass 
				<< setw(11) << energy << setw(18) << ndaughters << setw(16) << nparents << setw(10) << originalPDG << endl;
		}
		// get the information of original leptons
		LCObjectVec vecMCTL = navMCTL->getRelatedFromObjects(mcPart);
		FloatVec vecWgtMCTL = navMCTL->getRelatedFromWeights(mcPart);
		Int_t nMCTLRec = vecMCTL.size();
		Double_t chargeRec = 99;
		Double_t wgtRec = 0;
		if (vecMCTL.size() > 0) {
			ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(vecMCTL[0]);
			chargeRec = recPart->getCharge(); 
			wgtRec = vecWgtMCTL[0];
		}
		if (pdg == 11 && motherpdg == 0) {
			elepton1 = energy;
			cos1 = pv.CosTheta();
			phi1 = pv.Phi();
			lortzZ += lortz;
			lortzLep1 = lortz;
			chargeRec1 = chargeRec;
			wgtRec1 = wgtRec;
			nMCTL1 = nMCTLRec;
		}
		if (pdg == -11 && motherpdg == 0) {
			elepton2 = energy;
			cos2 = pv.CosTheta();
			phi2 = pv.Phi();
			lortzZ += lortz;
			lortzLep2 = lortz;
			chargeRec2 = chargeRec;
			wgtRec2 = wgtRec;
			nMCTL2 = nMCTLRec;
		}
		if (pdg == 13 && motherpdg == 0) {
			elepton3 = energy;
			cos3 = pv.CosTheta();
			phi3 = pv.Phi();
			lortzZ += lortz;
			lortzLep1 = lortz;
			chargeRec3 = chargeRec;
			wgtRec3 = wgtRec;
			nMCTL3 = nMCTLRec;
		}
		if (pdg == -13 && motherpdg == 0) {
			elepton4 = energy;
			cos4 = pv.CosTheta();
			phi4 = pv.Phi();
			lortzZ += lortz;
			lortzLep2 = lortz;
			chargeRec4 = chargeRec;
			wgtRec4 = wgtRec;
			nMCTL4 = nMCTLRec;
		}
		if (pdg == 25 && motherpdg == 0 && abs(daughterpdg) == 15) {
			nHtt +=1 ;
		}
		if (pdg == 25 && abs(daughterpdg) == 5) {
			nHbb +=1 ;
		}
		if (pdg == 25 && abs(daughterpdg) == 24) nHww++;
		if (abs(pdg) == 23 && abs(motherpdg) == 25 && (abs(daughterpdg) == 12 || abs(daughterpdg) == 14 || abs(daughterpdg) == 16)) nZvv++;
		if (abs(pdg) == 23 && abs(motherpdg) == 25 && abs(daughterpdg) < 10 && abs(daughterpdg) > 0) nZqq++;
		if (abs(pdg) == 23 && abs(motherpdg) == 25 && abs(daughterpdg) == 11) nZee++;
		if (abs(pdg) == 23 && abs(motherpdg) == 25 && abs(daughterpdg) == 13) nZmm++;
		if ((pdg == 11 || pdg == 13) && abs(motherpdg) == 23 && abs(mmotherpdg) == 25) {
			energyLep1MC = energy;
			cosLep1MC = pv.CosTheta();
		}
		if ((pdg == -11 || pdg == -13) && abs(motherpdg) == 23 && abs(mmotherpdg) == 25) {
			energyLep2MC = energy;
			cosLep2MC = pv.CosTheta();
		}
		if (pdg == 25 && abs(daughterpdg) == 21) nHgg++;
		if (pdg == 25 && abs(daughterpdg) == 23) nHzz++;
		if (pdg == 25 && abs(daughterpdg) == 4) nHcc++;
		Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
		if (nparents == 0 && pdg == 12 && ioverlay == 0) nn1++;
		if (nparents == 0 && pdg == 14 && ioverlay == 0) nn2++;
		if (nparents == 0 && pdg == 16 && ioverlay == 0) nn3++;
	}
	Double_t data_gen[15];
	data_gen[0] = nn1;
	data_gen[1] = nn2;
	data_gen[2] = nn3;
	data_gen[3] = nHbb;
	data_gen[4] = nHww;
	data_gen[5] = nZqq;
	data_gen[6] = nZvv;
	data_gen[7] = nZee;
	data_gen[8] = nZmm;
	data_gen[9] = nHgg;
	data_gen[10]= nHcc;
	data_gen[11]= nHzz;
	data_gen[12]= nHtt;
	hGen->Fill(data_gen);

	// -- Read out PFO information --
	LCCollection *colPFO = evt->getCollection(_colPFOs);
	if (!colPFO) {
		std::cerr << "No PFO Collection Found!" << std::endl;
		throw marlin::SkipEventException(this);
	}
	hStat->Fill(++selid);
	gCutName[(Int_t)selid] << "MCParticle and PandoraPFOs Collections found!" << ends;
	Int_t nPFOs = colPFO->getNumberOfElements();
	//  cerr << "Number of PFOs: " << nPFOs << endl;
	LCCollectionVec *pNewPFOsCollection    = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollectionVec *pLeptonsCollection    = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollectionVec *pIsoLeptonsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	pNewPFOsCollection->setSubset(true);
	pLeptonsCollection->setSubset(true);
	pIsoLeptonsCollection->setSubset(true);
	std::vector<lcio::ReconstructedParticle*> newPFOs;
	std::vector<lcio::ReconstructedParticle*> electrons;
	std::vector<lcio::ReconstructedParticle*> muons;
	std::vector<lcio::ReconstructedParticle*> isoleptons;

	static TNtupleD *hPfo = 0;
	if (!hPfo) {
		stringstream tupstr_pfo;
		tupstr_pfo << "ntracks:charge:mcpdg:motherpdg:deltae:mmotherpdg:ndaughters" << ":"
			<< "mcoriginal:energy:type:pid"                                  << ":"
			<< "totalcalenergy:momentum:ecalenergy:hcalenergy:coneenergy"    << ":"
			<< "elepton1:elepton2:elepton3:elepton4"                         << ":"
			<< "cos1:cos2:cos3:cos4:nmctl:mcwgt:ievt:irun"                   << ":"
			<< "nhits:ncones:nconechg:nconeneu:coneec:coneen:energylink"     << ":"
			<< "costheta:yokeenergy:energycor:momentumcor"                   << ":"
			<< "d0:z0:r0:deltad0:deltaz0:nsigd0:nsigz0:nsigr0:iov"              
			<< ends;
		hPfo = new TNtupleD("hPfo","",tupstr_pfo.str().data());
	}

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
	TLorentzVector pVis = TLorentzVector(0.,0.,0.,0.);
	for (Int_t i=0;i<nPFOs;i++) {
		ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
		LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(recPart);
		FloatVec vecWgtMCTL = navMCTL->getRelatedToWeights(recPart);
		Int_t mcpdg, motherpdg, mmotherpdg;
		Double_t mcwgt=0.;
		mcpdg      =  0;
		motherpdg  = -99999;
		mmotherpdg = -99999;
		Double_t deltaE = -99999.;
		Double_t energyLink = -99999.;
		Int_t mcoriginal   = 0;
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
			motherpdg = 0;
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
		Double_t data[100];
		data[0] = ntracks;
		data[1] = charge;
		data[2] = mcpdg;
		data[3] = motherpdg;
		data[4] = deltaE;
		data[5] = mmotherpdg;
		data[6] = mcndaughters;
		data[7] = mcoriginal;
		data[8] = energy;
		data[9] = itype;
		data[10]= pid;
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
				nHits += (*iCluster)->getCalorimeterHits().size();
				//nHits = calHits.size();
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
			//
			if (charge == -1 && abs(mcoriginal) == 11 && (mcpdg == +11 || mcpdg == +13) && motherpdg!=22 && energy > energyLep1) {
				energyLep1    = energy;
				energyLep1Fsr = energyCorr;
				cosLep1       = cosTheta;
				indexLep1     = i;
				momentumLep1  = momentum;
				nLep1Found++;
			}
			if (charge == +1 && abs(mcoriginal) == 11 && (mcpdg == -11 || mcpdg == -13) && motherpdg!=22 && energy > energyLep2) {
				energyLep2    = energy;
				energyLep2Fsr = energyCorr;
				cosLep2       = cosTheta;
				indexLep2     = i;
				momentumLep2  = momentum;
				nLep2Found++;
			}
			pVis += TLorentzVector(momentum,energy);
			// save the pfo information
			data[11] = totalCalEnergy;
			data[12] = momentumMagnitude;
			data[13] = ecalEnergy;
			data[14] = hcalEnergy;
			data[15] = coneEnergy;
			data[16] = elepton1;
			data[17] = elepton2;
			data[18] = elepton3;
			data[19] = elepton4;
			data[20] = cos1;
			data[21] = cos2;
			data[22] = cos3;
			data[23] = cos4;
			data[24] = nMCTL;
			data[25] = mcwgt;
			data[26] = _nEvt;
			data[27] = _nRun;
			data[28] = nHits;
			data[29] = nConePFOs;
			data[30] = nConeCharged;
			data[31] = nConeNeutral;
			data[32] = coneEC;
			data[33] = coneEN;
			data[34] = energyLink;
			data[35] = cosTheta;
			data[36] = yokeEnergy;
			data[37] = energyCorr;
			data[38] = momentumMagCorr;
			data[39] = d0;
			data[40] = z0;
			data[41] = r0;
			data[42] = deltad0;
			data[43] = deltaz0;
			data[44] = nsigd0;
			data[45] = nsigz0;
			data[46] = nsigr0;
			data[47] = iOverlay;
			if (iLepTune) {
				hPfo->Fill(data);
			}
			// select the leptons
			if (charge != 0 && 
					totalCalEnergy/momentumMagnitude       > fElectronCut1 && 
					totalCalEnergy/momentumMagnitude       < fElectronCut2 &&
					ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && 
					(momentumMagnitude > c0Electron + c1Electron*coneEC)) {
				if (nsigd0 > 50 || nsigz0 > 5) continue;   // contraint to primary vertex
				electrons.push_back(recPart);
			}
			if (charge != 0 && 
					totalCalEnergy/momentumMagnitude < fMuonCut1 && 
					yokeEnergy > fMuonCut3 && 
					(momentumMagnitude > c0Muon + c1Muon*coneEC)) {
				if (nsigd0 > 5 || nsigz0 > 5) continue;    // contraint to primary vertex
				muons.push_back(recPart);
			}
			// select the isolated leptons
			if (charge != 0 ){
			
			}	
		}
	}
	Int_t nelectrons = electrons.size();
	Int_t nmuons     = muons.size();
	//  cerr << "nelectrons: " << nelectrons << "  nmuons: " << nmuons << endl;
	Int_t    iLeptonType = 0;
	Double_t visEnergy = pVis.E();
	Double_t ptTot = pVis.Pt();
	Double_t pzTot = pVis.Pz();
	Double_t pxTot = pVis.Px();
	Double_t pyTot = pVis.Py();
	Double_t pTot  = pVis.P();

	//if (nelectrons < 2 && nmuons < 2) throw marlin::SkipEventException(this);
	hStat->Fill(++selid);
	gCutName[(Int_t)selid] << "nElectrons > 2 || nMuons > 2" << ends;

	//if (nelectrons+nmuons        > 2) throw marlin::SkipEventException(this);
	hStat->Fill(++selid);
	gCutName[(Int_t)selid] << "nElectrons + nMuons == 2" << ends;

	// find the lepton pair nearest to the Z mass
	std::vector<lcio::ReconstructedParticle*> electronZ;
	std::vector<lcio::ReconstructedParticle*> muonZ;
	std::vector<lcio::ReconstructedParticle*> leptonZ;
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

	Double_t energyLep1Rec=0.,energyLep2Rec=0.;
	Double_t energyLep1FSR=0.,energyLep2FSR=0.;
	Int_t mcOriginalLep1=0,mcOriginalLep2=0;
	if (muonZ.size() == 0 && electronZ.size() == 0) {
		//cerr << "no lepton pair candidate found!" << endl;
		//throw marlin::SkipEventException(this);
	}
	else if (muonZ.size() > 0) {
		iLeptonType = 2;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iMuon=muonZ.begin();iMuon<muonZ.end();iMuon++) {
			leptonZ.push_back(*iMuon);
			Int_t mcoriginal = 0;
			LCObjectVec vecMCTL = navMCTL->getRelatedToObjects((*iMuon));
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
	}
	else {
		iLeptonType = 1;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iElectron=electronZ.begin();iElectron<electronZ.end();iElectron++) {
			leptonZ.push_back(*iElectron);
			Int_t mcoriginal = 0;
			LCObjectVec vecMCTL = navMCTL->getRelatedToObjects((*iElectron));
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
	hStat->Fill(++selid);
	gCutName[(Int_t)selid] << "at least one lepton pair satisfies |Mass - MassZ| < 40" << ends;

	// create the ZPFOs collection
	LCCollectionVec           * pZCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	ReconstructedParticleImpl * recoZ        = new ReconstructedParticleImpl();
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
		Bool_t isLep=kFALSE;
		for (std::vector<lcio::ReconstructedParticle*>::const_iterator iLep=leptonZ.begin();iLep<leptonZ.end();++iLep) {
			if ((*iObj) == (*iLep)) isLep = kTRUE;
		}
		if (!isLep) pNewPFOsCollection->addElement(*iObj);
	}
	evt->addCollection(pNewPFOsCollection,    _colNewPFOs.c_str());
	evt->addCollection(pLeptonsCollection,    _colLeptons.c_str());
	evt->addCollection(pIsoLeptonsCollection, _colIsoLeptons.c_str());
	evt->addCollection(pZCollection,          _colZPFOs.c_str()  );

	Int_t nNewPFOs = pNewPFOsCollection->getNumberOfElements();
	// save the quantities to the hEvt
	Double_t vdata[100];
	vdata[ 0] = nMCP;
	vdata[ 1] = nPFOs;
	vdata[ 2] = nelectrons;
	vdata[ 3] = nmuons;
	vdata[ 4] = iLeptonType;
	vdata[ 5] = nHbb;
	vdata[ 6] = visEnergy;
	vdata[ 7] = pTot;
	vdata[ 8] = pxTot;
	vdata[ 9] = pyTot;
	vdata[10] = pzTot;
	vdata[11] = ptTot;
	vdata[12] = nHww;
	vdata[13] = nZqq;
	vdata[14] = nZee;
	vdata[15] = nZmm;
	vdata[16] = nn1;
	vdata[17] = nn2;
	vdata[18] = nn3;
	vdata[19] = nZvv;
	vdata[20] = nHgg;
	vdata[21] = nHcc;
	vdata[22] = nHzz;
	vdata[23] = nHtt;
	vdata[24] = massZ;
	vdata[25] = massZ_FSR;
	vdata[26] = energyLep1;
	vdata[27] = energyLep2;
	vdata[28] = energyLep1FSR;
	vdata[29] = energyLep2FSR;
	vdata[30] = energyLep1Rec;
	vdata[31] = energyLep2Rec;
	vdata[32] = energyLep1MC;
	vdata[33] = energyLep2MC;
	vdata[34] = cosLep1MC;
	vdata[35] = cosLep2MC;
	vdata[36] = mcOriginalLep1;
	vdata[37] = mcOriginalLep2;
	vdata[38] = massZ_Split;
	vdata[39] = nSplit;
	vdata[40] = nNewPFOs;
	hEvt->Fill(vdata);

	//-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

	last->cd();
}



void LeptonSelectionProcessor::check( LCEvent * evt ) { 
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void LeptonSelectionProcessor::end(){ 

	cerr << "LeptonSelectionProcessor::end()  " << name() 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< endl ;
	//
	cerr << "  =============" << endl;
	cerr << "   Cut Summary " << endl;
	cerr << "  =============" << endl;
	cerr << "   ll+4 Jet    " << endl;
	cerr << "  =============" << endl;
	cerr << endl
		<< "  -----------------------------------------------------------" << endl
		<< "   ID   No.Events    Cut Description                         " << endl
		<< "  -----------------------------------------------------------" << endl;

	for (int id=0; id<20 && gCutName[id].str().data()[0]; id++) {
		cerr << "  " << setw( 3) << id
			<< "  " << setw(10) << static_cast<int>(hStat->GetBinContent(id+1))
			<< "  : " << gCutName[id].str().data() << endl;
	}
	cerr << "  -----------------------------------------------------------" << endl;

}
