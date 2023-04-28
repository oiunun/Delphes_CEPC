#include "LGFastJetClustering.h"
#include "legendre.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "EVENT/LCFloatVec.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include <TMatrixDSymEigen.h> 
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include <UTIL/LCTOOLS.h>
#include "fastjet/ClusterSequence.hh"
//#include "fastjet/SISConePlugin.hh"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace Legendre ;

LGFastJetClustering aLGFastJetClustering ;

LGFastJetClustering::LGFastJetClustering() : Processor("LGFastJetClustering") {

	// Processor description
	_description = "FastJet Clustering ..." ;


	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputCollection",
			"Collection of reconstructed particles",
			_inputCollection,
			std::string("Unset") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputMCTruthMap",
			"Collection of reconstructed particles",
			_inputMCTruthMap,
			std::string("Unset") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"OutputCollection",
			"Name of collection with the found jets",
			_outputCollection,
			std::string("Unset") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"IsoLepCollection",
			"Name of collection with the Isolated Lepton",
			_outputIsoLepCol,
			std::string("IsoLepCol") );

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"RemainCollection",
			"Name of collection without the Isolated Lepton",
			_outputRemainCol,
			std::string("RemainPFOs") );

	registerProcessorParameter("Debug",
			"debug printout",
			_print,
			int(0)); 

	registerProcessorParameter("SaveNPZ",
			"Save NPZ file or not",
			_save_npz,
			int(0)); 

	registerProcessorParameter("Label",
			"event label",
			_label,
			int(5)); 

	registerProcessorParameter("Algorithm",
			"FastJet algorithm",
			sAlgorithm,
			std::string("antikt_algorithm")); 

	registerProcessorParameter("R",
			"R Parameter",
			_RPar,
			double(0.7)); 

	registerProcessorParameter("P",
			"P Parameter",
			_PPar,
			double(0.0)); 

	registerProcessorParameter("EjetMin",
			"Ejet",
			_eJet,
			double(10.0)); 

	registerProcessorParameter("PtMin",
			"PtMin Cut",
			_PtCut,
			double(0.0)); 

	registerProcessorParameter("NJets",
			"max nb of jets",
			_nJetMax,
			int(25)); 

	registerProcessorParameter("nPFOmin",
			"min nb of input PFOs",
			_nPFOmin,
			int(5)); 

	registerProcessorParameter("InclusiveExclusive",
			"Inclusive = 0; Exclusive != 0",
			_InclusiveExclusive,
			int(1)); 

	registerProcessorParameter("FillTree",
			"tuple",
			_fillTree,
			int(0)); 

	registerProcessorParameter("RemoveList",
			"partilces Removed before JetClustering",
			_RemoveList,
			vector<string>()); 

	registerProcessorParameter("EnergyCut",
			"Energy Threshold of particles",
			_EnergyCut,
			double(0)); 

}

void LGFastJetClustering::init() { 

	printParameters() ;

	Nv    =  5;
	Pm    = 50;
	Ne    =  0;
	dat.clear(); 

	if(_fillTree){

		if(_print>1)cout << "LGFastJetClustering: Making Tuples" << endl;

		// Declaration of Trees 
		_Etree = new TTree("ntp","jet information");
		_Etree->Branch("NRun",        &_nRun,       "NRun/I");
		_Etree->Branch("NEvt",        &_nEvt,       "NEvt/I");
		_Etree->Branch("Ecms",        &_eCMS,       "Ecms/F");
		_Etree->Branch("NJets",       &_nJets,      "NJets/I");
		_Etree->Branch("NJetsHE",     &_nJetsHE,    "NJetsHE/I");
		_Etree->Branch("RPar",        &_rp,         "RPar/D");

		_Ftree = new TTree("EFP","information saved for energy flow network");
		_Ftree->Branch("npar",      &_NPJ,          "npar/I");
		_Ftree->Branch("JetCharge", &_JetCharge,    "JetCharge/D");
		_Ftree->Branch("JetPDGID",  &_JetPDGID,     "JetPDGID/D");
		_Ftree->Branch("Energy",     _Energy,       "Energy[100]/D");
		_Ftree->Branch("Mom",        _Mom,          "Mom[100]/D");
		_Ftree->Branch("THETA",      _THETA,        "THETA[100]/D");
		_Ftree->Branch("PHI",        _PHI,          "PHI[100]/D");
		_Ftree->Branch("D0",         _D0,           "D0[100]/D");
		_Ftree->Branch("Z0",         _Z0,           "Z0[100]/D");
		_Ftree->Branch("PDGID",      _PDGID,        "PDGID[100]/I");
		_Ftree->Branch("Charge",     _Charge,       "Charge[100]/I");
		_Ftree->Branch("CosT",       _CosT,         "CosT[100]/D");
		_Ftree->Branch("PhiSinT",    _PhiSinT,      "PhiSinT[100]/D");
		_Ftree->Branch("logEn",      _logEn,        "logEn[100]/D");
		_Ftree->Branch("logEnF",     _logEnF,       "logEnF[100]/D");
		_Ftree->Branch("logMo",      _logMo,        "logMo[100]/D");
		_Ftree->Branch("logMoF",     _logMoF,       "logMoF[100]/D");
		_Ftree->Branch("delR",       _delR,         "delR[100]/D");

		_Ftree->Branch("lab_ll",    &_lab_ll,       "lab_ll/I");
		_Ftree->Branch("lab_cc",    &_lab_cc,       "lab_cc/I");
		_Ftree->Branch("lab_bb",    &_lab_bb,       "lab_bb/I");
		_Ftree->Branch("lab_uu",    &_lab_uu,       "lab_uu/I");
		_Ftree->Branch("lab_tt",    &_lab_tt,       "lab_tt/I");
		_Ftree->Branch("lab_gg",    &_lab_gg,       "lab_gg/I");
		_Ftree->Branch("lab_aa",    &_lab_aa,       "lab_aa/I");
		_Ftree->Branch("lab_ww",    &_lab_ww,       "lab_ww/I");
		_Ftree->Branch("lab_zz",    &_lab_zz,       "lab_zz/I");
		_Ftree->Branch("lab_az",    &_lab_az,       "lab_az/I");

		_Ctree = new TTree("CMB","jet information");
		_Ctree->Branch("npar",      &_NPJ,          "npar/I");
		_Ctree->Branch("TrueHmass", &_TrueHmass,    "TrueHmass/D");
		_Ctree->Branch("TrueZmass", &_TrueZmass,    "TrueZmass/D");
		_Ctree->Branch("Hmass",     &_Hmass,        "Hmass/D");
		_Ctree->Branch("Zmass",     &_Zmass,        "Zmass/D");
		_Ctree->Branch("llRec",     &_llRec,        "llRec/D");
		_Ctree->Branch("llMass",    &_llMass,       "llMass/D");
		_Ctree->Branch("hMass",     &_hMass,        "hMass/D");
		_Ctree->Branch("hRecMass",  &_hRecMass,     "hRecMass/D");
		_Ctree->Branch("BCL",       &_BCL,          "BCL/I");
		_Ctree->Branch("Energy",     _Energy,       "Energy[200]/D");
		_Ctree->Branch("logEn",      _logEn,        "logEn[200]/D");
		_Ctree->Branch("logEnF",     _logEnF,       "logEnF[200]/D");
		_Ctree->Branch("Mom",        _Mom,          "Mom[200]/D");
		_Ctree->Branch("logMo",      _logMo,        "logMo[200]/D");
		_Ctree->Branch("logMoF",     _logMoF,       "logMoF[200]/D");
		_Ctree->Branch("THETA",      _THETA,        "THETA[200]/D");
		_Ctree->Branch("PHI",        _PHI,          "PHI[200]/D");
		_Ctree->Branch("delR",       _delR,         "delR[200]/D");
		_Ctree->Branch("D0",         _D0,           "D0[200]/D");
		_Ctree->Branch("Z0",         _Z0,           "Z0[200]/D");
		_Ctree->Branch("PDGID",      _PDGID,        "PDGID[200]/I");
		_Ctree->Branch("Charge",     _Charge,       "Charge[200]/I");
		_Ctree->Branch("CosT",       _CosT,         "CosT[200]/D");
		_Ctree->Branch("SinT",       _SinT,         "SinT[200]/D");
		_Ctree->Branch("PhiSinT",    _PhiSinT,      "PhiSinT[200]/D");
		_Ctree->Branch("Sphericity",&_Sphericity,   "Sphericity/D");
		_Ctree->Branch("Aplanarity",&_Aplanarity,   "Aplanarity/D");
		_Ctree->Branch("C",         &_C,            "C/D");
		_Ctree->Branch("D",         &_D,            "D/D");
		_Ctree->Branch("Thrust",    &_Thrust,       "Thrust/D");
		_Ctree->Branch("Theta",     &_Theta,        "Theta/D");
		_Ctree->Branch("Phi",       &_Phi,          "Phi/D");
		_Ctree->Branch("Major",     &_Major,        "Major/D");
		_Ctree->Branch("Minor",     &_Minor,        "Minor/D");
		_Ctree->Branch("Thetap",    &_Thetap,       "Thrustp/D");
		_Ctree->Branch("Phip",      &_Phip,         "Phip/D");
		_Ctree->Branch("ThrEDMp",   &_ThrEDMp,      "ThrEDMp/D");
		_Ctree->Branch("FM",         _FM,           "FM[50]/D");
		_Ctree->Branch("FMI",        _FMI,          "FMI[50]/D");
		_Ctree->Branch("FMS",        _FMS,          "FMS[50]/D");
		_Ctree->Branch("Ci",         _Ci,           "Ci[28]/D");
		_Ctree->Branch("lab_ll",    &_lab_ll,       "lab_ll/I");
		_Ctree->Branch("lab_cc",    &_lab_cc,       "lab_cc/I");
		_Ctree->Branch("lab_bb",    &_lab_bb,       "lab_bb/I");
		_Ctree->Branch("lab_uu",    &_lab_uu,       "lab_uu/I");
		_Ctree->Branch("lab_tt",    &_lab_tt,       "lab_tt/I");
		_Ctree->Branch("lab_gg",    &_lab_gg,       "lab_gg/I");
		_Ctree->Branch("lab_aa",    &_lab_aa,       "lab_aa/I");
		_Ctree->Branch("lab_ww",    &_lab_ww,       "lab_ww/I");
		_Ctree->Branch("lab_zz",    &_lab_zz,       "lab_zz/I");
		_Ctree->Branch("lab_az",    &_lab_az,       "lab_az/I");


	}

	if(sAlgorithm=="kt_algorithm"){
		SetAlgorithm(fastjet::kt_algorithm);
	}else if(sAlgorithm=="cambridge_algorithm"){ 
		SetAlgorithm(fastjet::cambridge_algorithm);
	}else if(sAlgorithm=="antikt_algorithm"){ 
		SetAlgorithm(fastjet::antikt_algorithm);
	}else if(sAlgorithm=="ee_kt_algorithm"){ 
		SetAlgorithm(fastjet::ee_kt_algorithm);
	}else if(sAlgorithm=="genkt_algorithm"){ 
		SetAlgorithm(fastjet::genkt_algorithm);
	}else if(sAlgorithm=="ee_genkt_algorithm"){ 
		SetAlgorithm(fastjet::ee_genkt_algorithm);
	}else if(sAlgorithm=="plugin_algorithm"){ 
		SetAlgorithm(fastjet::plugin_algorithm);
	}else{
		printf("please chooose the proper algorithm, now it is %s\n", sAlgorithm.c_str());
		exit(1);
	}

	RemoveResonance = -9999; 

	if ( _RemoveList.size()>0){
		RemoveResonance = atof(_RemoveList[0].c_str());
		_Charged = 0 ;
		for( unsigned int i =1; i< _RemoveList.size();i++)
		{
			_PdgidList.push_back(atoi(_RemoveList[i].c_str()));
			if(_PdgidList[i-1]<0) _Charged=1;
		}
	}
}

void LGFastJetClustering::processRunHeader( LCRunHeader* run) { 

	_nRun = run->getRunNumber() ;

} 

void LGFastJetClustering::processEvent( LCEvent * evt ) { 


	_lab_ll= 0; 
	_lab_cc= 0; 
	_lab_bb= 0;  
	_lab_ee= 0; 
	_lab_uu= 0; 
	_lab_tt= 0; 
	_lab_gg= 0; 
	_lab_aa= 0; 
	_lab_zz= 0; 
	_lab_ww= 0; 
	_lab_az= 0; 
	_Hmass = 0.0; 
	_Zmass = 0.0; 
	_TrueHmass = 125.0; 
	_TrueZmass = 0.0; 

	int bcl=_label; 
	if( _label == 5 ) _lab_bb = 1; 
	if( _label == 4 ) _lab_cc = 1; 
	if( _label <= 3 ) _lab_ll = 1; 
	if( _label == 11) _lab_ee = 1; 
	if( _label == 13) _lab_uu = 1; 
	if( _label == 15) _lab_tt = 1; 
	if( _label == 21) _lab_gg = 1; 
	if( _label == 22) _lab_aa = 1; 
	if( _label == 23) _lab_zz = 1; 
	if( _label == 24) _lab_ww = 1; 
	if( _label == 45) _lab_az = 1; 


	Cluflag.setBit(LCIO::CHBIT_LONG);
	_jetsCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	_leptCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	_leftCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	_jetsCol->setSubset(Cluflag.getFlag()) ;
	_leptCol->setSubset(Cluflag.getFlag()) ;
	_leftCol->setSubset(Cluflag.getFlag()) ;

	_nRun = evt->getRunNumber();
	_nEvt = evt->getEventNumber();

	if(_print>0)cout <<"Run " << _nRun << " Evt " << _nEvt << endl;

	LCCollection* enflowcol  =evt->getCollection(_inputCollection);
	LCCollection* MCTruthMap =evt->getCollection(_inputMCTruthMap);
	LCRelationNavigator *navMCTL   = new LCRelationNavigator(MCTruthMap);

	int nenflow =  enflowcol->getNumberOfElements(); 
	if(_print>1)cout <<" # of tracks & clusters " << nenflow << endl;

	double px, py, pz, E;

	fastjet::JetAlgorithm algorithm = GetAlgorithm();
	if(_print>0)cout <<" Jet algorithm is  " << algorithm << endl;

	vector<fastjet::PseudoJet> input_particles;

	Double_t dMass=9999; 
	vector<int> plist; plist.clear();
	int pdgid1=0, pdgid2=0;
	int PDGID1=0, PDGID2=0;
	TLorentzVector p4p(0,0,0,0), p4m(0,0,0,0), p40(0,0,0,240);
	ReconstructedParticle *save_lep1=0, * save_lep2=0; 
	if ( RemoveResonance > 0 ) {
		if(_print>1)cout << " RemoveList " << _RemoveList[0] << " "<< _RemoveList[1] << " " << _RemoveList[2] << endl;
		for ( int i=0; i<nenflow-1 ; i++){
			ReconstructedParticle* enflow1 = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( i ));
			Double_t charge1 = enflow1->getCharge();
			if(_print>1)cout << " charge =  " << charge1 << endl;
			if( (int)fabs(charge1) != _Charged) continue;
			//
			LCObjectVec vecMCTL1            = navMCTL->getRelatedToObjects(enflow1);
			if (vecMCTL1.size() > 0) {
				for(unsigned int k=0; k< vecMCTL1.size(); k++){
					MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL1[k]);
					pdgid1=mcp->getPDG();
				}
				if(_print>1)cout << " pdgid1 " << pdgid1 << " "<< vecMCTL1.size() << endl;
			}
			int skip1=1;
			for(unsigned int s =0; s<_PdgidList.size(); s++)
				if( pdgid1 == _PdgidList[s]) skip1=0;
			if(skip1) continue;
			//
			Double_t energy1 = enflow1->getEnergy();
			TVector3 momentum1 = TVector3(enflow1->getMomentum());
			TLorentzVector p41 = TLorentzVector(momentum1,energy1);
			//
			for ( int j=i+1; j<nenflow ; j++){
				ReconstructedParticle* enflow2 = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( j ));
				Double_t charge2 = enflow2->getCharge();
				if( (int)fabs(charge2)!=_Charged) continue;
				if( fabs(charge1+charge2) > 1e-6) continue;
				//
				LCObjectVec vecMCTL2            = navMCTL->getRelatedToObjects(enflow2);
				if (vecMCTL2.size() > 0) {
					for(unsigned int k=0; k< vecMCTL2.size(); k++){
						MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL2[k]);
						pdgid2=mcp->getPDG();
					}
					if(_print>1)cout << " pdgid2 " << pdgid2 << " "<< vecMCTL2.size() << endl;
				}
				int skip2=1;
				for(unsigned int s =0; s<_PdgidList.size(); s++)
					if( pdgid2 == _PdgidList[s]) skip2=0;
				if(skip2) continue;
				//
				Double_t energy2 = enflow2->getEnergy();
				TVector3 momentum2 = TVector3(enflow2->getMomentum());
				TLorentzVector p42 = TLorentzVector(momentum2,energy2);
				double dmass = fabs( (p41+p42).M() - RemoveResonance ); 
				if ( dmass < dMass){
					save_lep1= enflow1;
					save_lep2= enflow2;
					p4p = p41, p4m = p42; 
					_llMass = ( p41 + p42 ).M() ;
					_llRec  = ( p40 - p41 - p42 ).M(); 
					dMass=dmass;
					PDGID1=pdgid1;
					PDGID2=pdgid2;
					plist.clear();
					plist.push_back(i);
					plist.push_back(j);
				}
			}
		}
	}

	if(save_lep1!=0) _leptCol->addElement(save_lep1);
	if(save_lep2!=0) _leptCol->addElement(save_lep2);
	for ( int i=0; i<nenflow ; i++){
		ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( i ));
		if(save_lep1!=enflow && save_lep2!=enflow) _leftCol->addElement(enflow);
	}
	//
	if(_print>1)cout << "PDGID1 & PDGID2 " << PDGID1<<" "<<PDGID2 <<" dMass = "<<dMass<< endl;
	if(_print>1)cout << "# of RemoveList is " << plist.size() << endl;
	//
	//
	//
	vector<ReconstructedParticle*> cmb4; cmb4.clear(); 
	TLorentzVector v4mc(0,0,0,0);
	TLorentzVector v4sm(0,0,0,0);
	for ( int ienflow=0; ienflow<nenflow ; ienflow++){
		ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( ienflow ));
		MCParticle* mcp = NULL; 
		if( _EnergyCut <  enflow->getEnergy() ) {
			//
			LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(enflow);
			int pdgid=0; 
			if (vecMCTL.size() > 0) {
				mcp = dynamic_cast<MCParticle *>(vecMCTL[0]);
				if(mcp)pdgid=abs(mcp->getPDG());
			}
			if ( (pdgid ==11 || pdgid==13 || pdgid== 211 || pdgid==321 || pdgid == 2212) && enflow->getEnergy()<_EnergyCut ) continue;
			else if ( (pdgid ==22)  && enflow->getEnergy()<0.5 ) continue;
			else if ( (pdgid ==130 || pdgid == 2112) && enflow->getEnergy()<1.0 ) continue;
			//
			int skip = 0;
			for( unsigned int i = 0; i < plist.size(); i++){
				if ( ienflow == plist[i] ) skip++;	
			}
			if( skip == 0 ){
				px = enflow->getMomentum()[0];
				py = enflow->getMomentum()[1];
				pz = enflow->getMomentum()[2];
				E  = enflow->getEnergy();
				TLorentzVector v4(mcp->getMomentum(), mcp->getEnergy());
				v4mc += v4;
				v4sm += TLorentzVector(px,py,pz,E);

				fastjet::PseudoJet thisPtc(px,py,pz,E);
				thisPtc.set_user_index(ienflow);
				if(_print>2) printf("id, px, py, pz, E= %5d, %10.4f, %10.4f, %10.4f, %10.4f\n", enflow->id(),px, py, pz, E);

				input_particles.push_back(thisPtc);
				cmb4.push_back(enflow);
			}
		}
	}

	_TrueZmass = v4mc.M();
	_Zmass = v4sm.M();
	if(_print>0) printf ("Z mass is %8.2f %8.2f\n",_TrueZmass, _Zmass); 

	if(_print>0) printf("%3d PFOs used for jet-clustering\n", (int)input_particles.size());
	double _ymin[8];
	memset(_ymin, 0, sizeof(double)*8 );

	vector<fastjet::PseudoJet> sortedInputs = fastjet::sorted_by_E(input_particles);
	int inp = input_particles.size(), out=0;	
	vector<ReconstructedParticleImpl*> trash;
	float momentum[3], energy, mass;
	if( (int)input_particles.size() > _nPFOmin ) {
		_rp = _RPar;
		_pp = _PPar;

		fastjet::Strategy strategy = fastjet::Best;

		fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;


		fastjet::JetDefinition* jet_def =NULL;

		if(algorithm==fastjet::kt_algorithm || algorithm==fastjet::antikt_algorithm || algorithm==fastjet::cambridge_algorithm){
			jet_def = new fastjet::JetDefinition(algorithm, _rp, recomb_scheme, strategy); 
		}else if(algorithm==fastjet::cambridge_algorithm){
			jet_def = new fastjet::JetDefinition(algorithm, _rp, recomb_scheme, strategy);
		}else if(algorithm==fastjet::ee_kt_algorithm){ 
			jet_def = new fastjet::JetDefinition(algorithm); 
		}else if(algorithm==fastjet::genkt_algorithm || algorithm==fastjet::ee_genkt_algorithm ){ 
			jet_def = new fastjet::JetDefinition(algorithm, _rp, _pp, recomb_scheme, strategy); 
		}else if (algorithm==fastjet::plugin_algorithm ){
			double cone_radius = _rp;
			double overlap_threshold = _pp;
			//fastjet::SISConePlugin siscone(cone_radius, overlap_threshold);
			//jet_def = new fastjet::JetDefinition(& siscone);
		}

		fastjet::ClusterSequence cs(input_particles, *jet_def);

		vector<fastjet::PseudoJet> JetsVec ;

		if( _InclusiveExclusive ) JetsVec = cs.exclusive_jets(_nJetMax);
		else                      JetsVec = cs.inclusive_jets(_PtCut  );

		vector<fastjet::PseudoJet> sortedJets = fastjet::sorted_by_E(JetsVec);


		for(int iy=1; iy<9;iy++){	
			_ymin[iy-1] = cs.exclusive_ymerge (iy); // the ymin corresponding to the recombination that went from iy+1 to iy jets
			if(_print>1)printf("%3d,  ymin = %12.5f\n", iy, _ymin[iy-1] );
		}	


		int nmx = sortedJets.size();
		if( _InclusiveExclusive>0 && nmx > _nJetMax) nmx = _nJetMax;

		_nJets = sortedJets.size();

		delete jet_def; 

		//analyze( sortedJets) ;

		if(_print>0)printf("LGFastJetClustering: Nb of Jets %3d\n", (int)sortedJets.size() );
		_nJetsHE=0;
		for(int ij=0; ij<_nJets;ij++){

			vector<ReconstructedParticle*> cnn4; cnn4.clear(); 
			momentum[0] = sortedJets[ij].px();
			momentum[1] = sortedJets[ij].py();
			momentum[2] = sortedJets[ij].pz();
			energy      = sortedJets[ij].e();
			mass        = sortedJets[ij].m();
			TLorentzVector pJet(sortedJets[ij].px(), sortedJets[ij].py(), sortedJets[ij].pz(), sortedJets[ij].e());


			vector<double> npzEVT(Nv*Pm);
			if( _save_npz > 0 ){
				for(unsigned int i=0; i<Nv*Pm ; i++){
					npzEVT[i] = 0.0; 
					dat.push_back(0.0); 
				}
			}

			if(_print>1) printf( "Jet %2d  Energy = %10.4f\n", ij,  energy );

			if(energy>_eJet && ij < 10 ){

				vector<fastjet::PseudoJet> jetConstituents = cs.constituents(sortedJets[ij]);
				out+= jetConstituents.size();

				_nJetsHE++;

				map<MCParticle *, int > parentmap;
				if( _InclusiveExclusive==0 || ij<_nJetMax){
					ReconstructedParticleImpl* Jet = new ReconstructedParticleImpl;
					_NPJ = jetConstituents.size(); 
					for(unsigned ip = 0; ip < jetConstituents.size(); ip++){
						if(jetConstituents[ip].user_index()>=0 && jetConstituents[ip].user_index() < enflowcol->getNumberOfElements()){
							ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt(jetConstituents[ip].user_index()));
							cnn4.push_back(enflow);
							//
							LCObjectVec vecMCTL1  = navMCTL->getRelatedToObjects(enflow);
							MCParticle * mcparent = NULL, *mcp = NULL;
							if (vecMCTL1.size() > 0) {
								for(unsigned int k=0; k< vecMCTL1.size(); k++){
									mcp = dynamic_cast<MCParticle *>(vecMCTL1[k]);
									mcparent = FindParton(mcp);
								}
							}
							if( mcparent != NULL ){
								map<MCParticle *, int>::iterator it;
								it = parentmap.find(mcparent);
								if( it != parentmap.end() ) 
									parentmap[mcparent] += 1;
								else
									parentmap[mcparent]  = 1;
							}
							//
							Jet->addParticle(enflow);
							if ( enflow->  getTracks().size()>0){
								if(fabs(enflow->getCharge())>0.1) Jet->addTrack  (enflow->  getTracks()[0]);
							}
							if ( enflow->  getClusters().size()>0){
								if(fabs(enflow->getCharge())<0.1) Jet->addCluster(enflow->getClusters()[0]);
							}
							if(_print>2) cout << "add to particle list  " << enflow->id() << endl;
						}
					}

					vector<ReconstructedParticle*> sorted_cnn= sorted_by_E(cnn4);
					double OPhi=-99,OTheta=-99;
					unsigned int Np = 0; 	
					for(int ik=0; ik<100;ik++){
						_Energy  [ik]  =   0.;
						_Mom     [ik]  =   0.; 
						_logEn   [ik]  = -999; 
						_logEnF  [ik]  = -999; 
						_logMo   [ik]  = -999; 
						_logMoF  [ik]  = -999; 
						_delR    [ik]  = -999; 
						_THETA   [ik]  = -999; 
						_PHI     [ik]  = -999; 
						_D0      [ik]  =    0; 
						_Z0      [ik]  =    0; 
						_CosT    [ik]  = -999; 
						_PhiSinT [ik]  = -999; 
						_Charge  [ik]  = -999; 
						_PDGID   [ik]  =    0; 
					}
					for ( unsigned int ip =0; ip < sorted_cnn.size(); ip++)
					{
						if( ip>100) break;
						LCObjectVec vecMCTL1  = navMCTL->getRelatedToObjects(sorted_cnn[ip]);
						MCParticle *mcp = NULL;
						if (vecMCTL1.size() > 0) {
							mcp = dynamic_cast<MCParticle *>(vecMCTL1[0]);
						}
						TLorentzVector p40(sorted_cnn[ip]->getMomentum(),sorted_cnn[ip]->getEnergy());
						if( 0 == ip ){
							//double * angle2 = angle(p40); // leading particle
							double * angle2 = angle(pJet);   // jet axis 
							OTheta  = *(angle2+0);
							OPhi    = *(angle2+1);
						}
						TLorentzVector p41 = erout4(p40, OPhi,OTheta);
						_Energy [ip]  = p41.T();
						_logEn  [ip]  = std::log(p41.T());
						_logEnF [ip]  = std::log(p41.T()/pJet.T());
						_Mom    [ip]  = p41.Rho(); 
						_logMo  [ip]  = std::log(p41.Rho());
						_logMoF [ip]  = std::log(p41.Rho()/pJet.Rho());
						_Charge [ip]  = mcp->getCharge(); 
						_CosT   [ip]  = p41.CosTheta(); 
						_THETA  [ip]  = p41.Theta(); 
						_PHI    [ip]  = p41.Phi(); 
						_delR   [ip]  = _THETA[ip] * _THETA[ip] +_PHI[ip] * _PHI[ip]; 
						_delR   [ip]  = sqrt(_delR[ip]);  
						_PhiSinT[ip]  = p41.Phi()*sin(p41.Theta());
						_PDGID  [ip]  = mcp->getPDG(); 

						double d0  = 0, z0 = 0; 
						if ( mcp && (abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13 || abs(mcp->getPDG())==211 || abs(mcp->getPDG())==321 ||  abs(mcp->getPDG())==2212) ){
							const double *v = mcp->getVertex(); 
							d0 = pow(v[0]*v[0] +v[1]*v[1],0.5);  
							z0 = v[2];  
							if( ( d0*d0+z0*z0 ) < 10000 ){
								_D0    [ip]  = d0; 
								_Z0    [ip]  = z0; 
							}
						}
						//
						if( ip < Pm && _save_npz>0 ){
							npzEVT[ Np*Nv + 0 ] =  p41.T();
							npzEVT[ Np*Nv + 1 ] =  p41.Theta();
							npzEVT[ Np*Nv + 2 ] =  p41.Phi();
							npzEVT[ Np*Nv + 3 ] =  mcp->getPDG() ;
							npzEVT[ Np*Nv + 4 ] =  d0      ;
						}
						Np++; 
					}

					Jet->setMomentum(momentum);
					Jet->setEnergy(energy);
					Jet->setMass(mass);
					Jet->setType(4);
					_jetsCol->addElement(Jet);
					trash.push_back(Jet);	
					//
					map<MCParticle *, int>::iterator it;
					int ic=-99;
					for(it = parentmap.begin();it != parentmap.end();it++){
						if ( it->second > ic && it->first->getPDG()!=0  ){
							bcl   = 5-abs(it->first->getPDG());
							ic    = it->second;	
							_JetCharge = it->first->getCharge()*3.0+3.0;
							_JetPDGID  = it->first->getPDG();
						}
					}
					parentmap.clear();
					bcl = _label; 
					if(_save_npz >0 ) {
						if( bcl>0 ){
							Ne++;
							tag.push_back(bcl); 
							for(unsigned int i=0; i<Nv*Pm ; i++){
								dat[Nv*Pm*(Ne-1)+i]=npzEVT[i]; 
							}
						}
						npzEVT.clear(); 
					}
				}
			}
			if(_fillTree  ){
				_Ftree->Fill();
			}
		}

		if(_print>0) printf("input  %3d  output = %3d\n", inp, out);
		if(_fillTree){
			_Etree->Fill();
		}    
	}

	if( _fillTree ){
		//
		_BCL = bcl; 
		TMatrixDSym sphMat(3);
		sphMat(0,0) = 0;
		sphMat(0,1) = 0;
		sphMat(0,2) = 0;
		sphMat(1,0) = 0;
		sphMat(1,1) = 0;
		sphMat(1,2) = 0;
		sphMat(2,0) = 0;
		sphMat(2,1) = 0;
		sphMat(2,2) = 0;
		//
		double OPhi=-99,OTheta=-99;
		for(int ij=0; ij<200;ij++){
			_Energy [ij]  =  0.;
			_Mom    [ij]  =  0.; 
			_logEn  [ij]  = -999; 
			_logEnF [ij]  = -999; 
			_logMo  [ij]  = -999; 
			_logMoF [ij]  = -999; 
			_delR   [ij]  = -999; 
			_PHI    [ij]  = -999; 
			_THETA  [ij]  = -999; 
			_PDGID  [ij]  =  0; 
		   _Charge [ij]  = -999; 
			_D0     [ij]  =  0; 
			_Z0     [ij]  =  0; 
			_PhiSinT[ij]  = -999;
			_CosT   [ij]  = -999; 
			_SinT   [ij]  = -999; 
		}
		vector<ReconstructedParticle*> sorted_cmb= sorted_by_E(cmb4);
		TLorentzVector pevt(0,0,0,0);
		for(int ij=0; ij<inp;ij++){
			pevt += TLorentzVector(sorted_cmb[ij]->getMomentum(),sorted_cmb[ij]->getEnergy());
		}
		_hMass    = pevt.M();
		_hRecMass = ( p40 - pevt ).M();
		_NPJ = sorted_cmb.size(); 
		for(int ij=0; ij<inp;ij++){
			TLorentzVector p40(sorted_cmb[ij]->getMomentum(),sorted_cmb[ij]->getEnergy());
			//
			sphMat(0,0) += p40.X()*p40.X();
			sphMat(0,1) += p40.X()*p40.Y();
			sphMat(0,2) += p40.X()*p40.Z();
			sphMat(1,0) += p40.Y()*p40.X();
			sphMat(1,1) += p40.Y()*p40.Y();
			sphMat(1,2) += p40.Y()*p40.Z();
			sphMat(2,0) += p40.Z()*p40.X();
			sphMat(2,1) += p40.Z()*p40.Y();
			sphMat(2,2) += p40.Z()*p40.Z();
			if( ij>200) break;

			LCObjectVec vecMCTL1  = navMCTL->getRelatedToObjects(sorted_cmb[ij]);
			MCParticle *mcp = NULL;
			if (vecMCTL1.size() > 0) {
				mcp = dynamic_cast<MCParticle *>(vecMCTL1[0]);
			}
			if( 0 == ij ){
				double * angle2 = angle(pevt);
				OTheta  = *(angle2+0);
				OPhi    = *(angle2+1);
			}
			TLorentzVector p41 = erout4(p40, OPhi, OTheta);
			_THETA[ij]   = p41.Theta(); 
			_PHI[ij]     = p41.Phi();
			_Energy[ij]  = p41.E();
			_Mom   [ij]  = p41.Rho();
			_CosT[ij]    = cos(p41.Theta());
			_SinT[ij]    = sin(p41.Theta());
			_PhiSinT[ij] = p41.Phi()*sin(p41.Theta());
			//
			_logEn  [ij]  = std::log(p41.T());
			_logEnF [ij]  = std::log(p41.T()/pevt.T());
			_logMo  [ij]  = std::log(p41.Rho());
			_logMoF [ij]  = std::log(p41.Rho()/pevt.Rho());
			_delR   [ij]  = _THETA[ij] * _THETA[ij] +_PHI[ij] * _PHI[ij]; 
			_delR   [ij]  = sqrt(_delR[ij]);  
			//
			if( mcp ){
				_PDGID [ij]  = mcp->getPDG();
				_Charge[ij]  = mcp->getCharge(); 
			}
			TrackVec   trkvec = sorted_cmb[ij]->getTracks();
			if( trkvec.size()>0){
				if ( mcp && (abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13 || abs(mcp->getPDG())==211 || abs(mcp->getPDG())==321 ||  abs(mcp->getPDG())==2212 ) ) {
					double d0 = trkvec[0]->getD0();
					double z0 = trkvec[0]->getZ0();
					if( ( d0*d0+z0*z0 ) < 10000 ){
						_D0[ij] = d0;
						_Z0[ij] = z0;
					}
				}
				//if( 0 >= ij ) printf("%4d, %8.3f %8.3f %8.3f %8.3f, %8.0f, %9.3f %9.3f\n", ij, p40.Rho(), p41.Rho(), p41.Phi(), p41.Theta(), _PDGID[ij], _D0[ij], _Z0[ij]);
			}
		}

		double norm = sphMat(0,0)+sphMat(1,1)+sphMat(2,2);
		sphMat *= 1.5/norm;
		//
		TMatrixDSymEigen EigMat(sphMat);
		TVectorD eig   = EigMat.GetEigenValues();
		_Sphericity =  eig(1)+eig(2);
		_Aplanarity =  eig(2);
		_C          =  (eig(0)*eig(1)+eig(0)*eig(2)+eig(1)*eig(2))*2;
		_D          =  (eig(0)*eig(1)*eig(2))*2*9;
		//
		thrust(sortedInputs);
		fox_wolfram(sortedInputs);
		//if( bcl <2) 
		_Ctree->Fill();
	}

	_jetsCol->parameters().setValue( "RPar" , (float)_rp );
	_jetsCol->parameters().setValue( "PPar" , (float)_pp );
	_jetsCol->parameters().setValue( "NJets", (float)_nJets );
	for( int iy=1; iy<9; iy++){
		stringstream yname, dname;
		yname << "y"<<iy<<iy+1;
		dname << "d"<<iy<<iy+1;
		_jetsCol->parameters().setValue( yname.str(), (float)_ymin[iy-1] );
	}

	_jetsCol->setDefault   ( true   ) ; // only true works, false has track/cluster but without particle
	_jetsCol->setSubset    ( false  ) ; // can make reasonble slcio file with PFO collection/LCRelation !!! 
	_jetsCol->setTransient ( false  ) ;

	evt->addCollection(_jetsCol ,_outputCollection) ;
	evt->addCollection(_leptCol ,_outputIsoLepCol) ;
	evt->addCollection(_leftCol ,_outputRemainCol) ;

	//if(_fillTree) processJets(evt, _jetsCol);
	//FreeDelAll(trash);	
	delete navMCTL;
}

void  LGFastJetClustering::check( LCEvent * evt ) { 
}

void  LGFastJetClustering:: processJets( LCEvent * evt , LCCollectionVec* jets ) { 
}
void LGFastJetClustering:: fox_wolfram(const vector<fastjet::PseudoJet> sortedInputs ) 
{
	cout.precision(5) ;
	vector<TLorentzVector> p4list4FM;
	const unsigned inp = sortedInputs.size();
	double sqrts=0; 
	for( unsigned n =0;n<50; n++){
		_FM [n] = 0.0; 
		_FMI[n] = 0.0; 
		_FMS[n] = 0.0; 
	}
	if(  sortedInputs.size()>1){
		p4list4FM.clear();
		for (unsigned int i = 0; i<inp; i++) {
			TLorentzVector v4(sortedInputs[i].px(),sortedInputs[i].py(), sortedInputs[i].pz(), sortedInputs[i].e());
			sqrts+=sortedInputs[i].e();
			p4list4FM.push_back(v4);
		}	
		sqrts = sqrts*sqrts; 
		for (unsigned int i = 0; i<inp; i++) {
			double eni = p4list4FM[i].T(); 
			double csi = cos(p4list4FM[i].Theta()); 
			double ssi = sin(p4list4FM[i].Theta()); 
			double phi = p4list4FM[i].Phi(); 
			for (unsigned int j = 0; j<inp; j++) {
				double enj = p4list4FM[j].T(); 
				double csj = cos(p4list4FM[j].Theta()); 
				double ssj = sin(p4list4FM[j].Theta()); 
				double phj = p4list4FM[j].Phi(); 
				double x = csi*csj+ssi*ssj*cos(phi-phj); 
				for( unsigned n =0; n<50; n++){
					double pn = Pn(n,x);
					//printf("P%2d(%7.4f) = %8.4f\n",n,x,pn);
					pn = pn*(eni*enj);
					if(i==j) _FMS[n] += pn ; 
					else     _FMI[n] += pn; 
					_FM[n] += pn ; 
				}
			}
		}
		for( unsigned n =0;n<50; n++){
			_FM [n] /= sqrts; 
			_FMI[n] /= sqrts; 
			_FMS[n] /= sqrts; 
		}
	}

}	
//**********************************************
vector<TVector3> p3list4thrust;
void LGFastJetClustering::thrust( const vector<fastjet::PseudoJet> sortedInputs ) {

	if(  sortedInputs.size()>1){
		p3list4thrust.clear();
		TVector3 p3(0,0,0);
		for (unsigned int i = 0; i<sortedInputs.size(); i++) {
			TLorentzVector trkVec(sortedInputs[i].px(),sortedInputs[i].py(), sortedInputs[i].pz(),sortedInputs[i].e());
			TVector3 v3(trkVec.Vect());
			p3list4thrust.push_back(v3);
			p3 += v3;
		}	
		//
		void FCN0A(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
		void FCNPA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
		TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
		gMinuit->SetFCN(FCN0A);
		//
		Int_t    ierflg = 0;
		Double_t arglist[10], val[10], err[10];
		//
		arglist[0] =-1.0;
		gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
		arglist[0] = 1.0;
		gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
		arglist[0] = 2.0;
		gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
		//
		// Set starting values and step sizes for parameters
		//
		Double_t vstart[2] = {p3.Theta(), p3.Phi() };
		Double_t   step[2] = {    0.0100,  0.01000 };
		gMinuit->mnparm(0, "Theta",  vstart[0], step[0], 0.000, 0.000, ierflg);
		gMinuit->mnparm(1, "Phi",    vstart[1], step[1], 0.000, 0.000, ierflg);
		//
		// Now ready for minimization step
		//
		arglist[0] = 5000;
		arglist[1] = 1.e-8;
		gMinuit->mnexcm("SIMPLEX" , arglist , 2, ierflg); // mimimization here ... 
		//Print results
		Double_t amin,edm,errdef;
		Int_t nvpar,nparx,icstat;
		gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		_ThrEDM =  edm;
		_Thrust =-amin;
		//
		//
		TString chnam;
		Double_t xlolim, xuplim;
		Int_t    iuint;
		for(Int_t p = 0; p < 2; p++)
		{ 
			gMinuit->mnpout(p, chnam, val[p], err[p], xlolim, xuplim, iuint);
		}
		val[0] = fmod(val[0]+CLHEP::twopi+CLHEP::twopi, CLHEP::pi);
		val[1] = fmod(val[1]+CLHEP::twopi+CLHEP::twopi+CLHEP::pi, CLHEP::twopi) - CLHEP::pi;
		_Theta = val[0];
		_Phi   = val[1];
		TVector3 n0(sin(_Theta)*cos(_Phi), sin(_Theta)*sin(_Phi), cos(_Theta));
		//printf("theta = %10.5f --> %10.5f; phi = %10.5f --> %10.5f\n", vstart[0], val[0], vstart[1], val[1]);
		delete gMinuit; gMinuit=NULL;
		//
		TVector3 pt(0,0,0);
		for (unsigned int i = 0; i<p3list4thrust.size(); i++) {
			pt += p3list4thrust[i];
		}	
		gMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 4 params
		gMinuit->SetFCN(FCNPA);
		//
		arglist[0] =-1.0;
		gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);
		arglist[0] = 1.0;
		gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
		arglist[0] = 2.0;
		gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
		//
		// Set starting values and step sizes for parameters
		//
		TVector3 nz(1,1,1);
		TVector3 ng=n0.Cross(pt.Unit());
		Double_t vstartp[4] = {  ng.Theta(), ng.Phi(), val[0], val[1] };
		Double_t   stepp[4] = {      0.1000,   0.1000, 0.0000, 0.0000 };
		gMinuit->mnparm(0, "Thetap",  vstartp[0], stepp[0], 0.000, 0.000, ierflg);
		gMinuit->mnparm(1, "Phip",    vstartp[1], stepp[1], 0.000, 0.000, ierflg);
		gMinuit->mnparm(2, "Theta",   vstartp[2], stepp[2], 0.000, 0.000, ierflg);
		gMinuit->mnparm(3, "Phi",     vstartp[3], stepp[3], 0.000, 0.000, ierflg);
		//
		// Now ready for minimization step
		//
		arglist[0] = 5000;
		arglist[1] = 1.e-8;
		gMinuit->mnexcm("SIMPLEX", arglist , 2, ierflg); // mimimization here ... 
		//gMinuit->mnexcm("MINI"   , arglist , 2, ierflg); // mimimization here ... 
		//Print results
		gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		_ThrEDMp =   edm;
		_Major   = -amin;
		//
		for(Int_t p = 0; p < 2; p++)
		{ 
			gMinuit->mnpout(p, chnam, val[2+p], err[p], xlolim, xuplim, iuint);
		}
		//
		val[2] = fmod(val[2]+CLHEP::twopi+CLHEP::twopi, CLHEP::pi);
		val[3] = fmod(val[3]+CLHEP::twopi+CLHEP::twopi+CLHEP::pi, CLHEP::twopi) - CLHEP::pi;
		_Thetap = val[2];
		_Phip   = val[3];
		TVector3 n1(sin(_Thetap)*cos(_Phip),sin(_Thetap)*sin(_Phip), cos(_Thetap));
		//	
		double minor = 0, sump=0;
		TVector3 nm=n1.Cross(n0);
		for (unsigned int i = 0; i<p3list4thrust.size(); i++) {
			TVector3 v = p3list4thrust[i];
			double  p = v.Mag();
			minor    += fabs(v.Dot(nm));
			sump     += p;	
		}
		_Minor   = minor/sump;
		//
		//printf("thetap= %10.5f --> %10.5f; phip= %10.5f --> %10.5f: %11.8f\n", ng.Theta(), val[2], ng.Phi(), val[3], n0.Dot(n1));
		delete gMinuit; gMinuit=NULL;
		//
	}
}
void FCN0A(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const double theta = par[0];
	const double phi   = par[1];
	TVector3 n0(sin(theta)*cos(phi),sin(theta)*sin(phi), cos(theta));
	//calculate T 
	double chisq = 0, sump=0;
	for (unsigned int i = 0; i<p3list4thrust.size(); i++) {
		TVector3 v = p3list4thrust[i];
		double  p = v.Mag();
		chisq    += fabs(v.Dot(n0));
		sump     += p;	
	}
	//
	f = -chisq/sump;
}
//**********************************************
//**********************************************
void FCNPA(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const double theta = par[0];
	const double phi   = par[1];
	const double thetap= par[2];
	const double phip  = par[3];
	TVector3 n0(sin(theta )*cos(phi ),sin(theta )*sin(phi ), cos(theta ));
	TVector3 n1(sin(thetap)*cos(phip),sin(thetap)*sin(phip), cos(thetap));
	double penalty = n0.Dot(n1);
	//calculate T 
	double chisq = 0, sump=0;
	for (unsigned int i = 0; i<p3list4thrust.size(); i++) {
		TVector3 v = p3list4thrust[i];
		double  p = v.Mag();
		chisq    += fabs(v.Dot(n0));
		sump     += p;	
	}
	//
	f = -chisq/sump + penalty*penalty*999;
}


void LGFastJetClustering::end(){ 
	std::cout << "LGFastJetClustering::end()  "  << " processed " << _nEvt << " events in " << _nRun << " runs " << std::endl ;
	//
	if ( _save_npz > 0 ) {
		cnpy::npz_save("out.npz","X",&dat[0],{Ne,  Pm,  Nv}, "w"); 
		cnpy::npz_save("out.npz","Y",&tag[0],{Ne          }, "a"); 
	}
	dat.clear(); 
	tag.clear();
	//
	if(_fillTree){
		if(_print>1)cout << "LGFastJetClustering: Saving Tuples" << endl;
	}
}
//
//**********************************************
//**********************************************
MCParticle * LGFastJetClustering::FindParton(MCParticle *mcp)
{
	MCParticle *ret = NULL;

	if ( mcp->getParents().size() == 0 )      // original particle, no parent  
	{
		if ( abs(mcp->getPDG())<=6 ){
			return mcp;
		}
		else return NULL;
	} 
	else if ( mcp->getParents().size() == 1 ) // Only one parent  
	{
		MCParticle *mcpa = mcp->getParents()[0];
		//
		if ( mcp->getPDG()==mcpa->getPDG() )
		{  
			// mother is same particle
			ret = FindParton(mcpa);	
		}
		else if ( abs(mcp->getPDG())<=6 && ( abs(mcpa->getPDG())>=23 && abs(mcpa->getPDG())<=25 ) )
		{
			// parent is W/Z and itself is a quark 
			ret = mcp;
		}
		else if ( mcp->getPDG()==21  && mcpa->getPDG()==25 )
		{
			// parent is Higgs and itself is a gluon 
			ret = mcp;
		}
		else if ( abs(mcp->getPDG())<= 6 && mcpa->getPDG()==21 )
		{
			// parent is gluon and itself is a quark, this kind of object will be neglected   
			ret = NULL;
		}	
		else if ( abs(mcpa->getPDG())<= 6 && mcp->getPDG()==21 )
		{
			// parent is quark and itself is a gluon 
			ret = FindParton(mcpa);	
		}	
		else if ( abs(mcpa->getPDG())>= 81 && abs(mcpa->getPDG())<=100 && mcpa->getParents().size()>=2 ) 
		{
			// parent is intermediate particle; need to check the MANY(>1) grandparents
			if ( abs(mcp->getPDG()) <= 6 || ( abs(mcp->getPDG()) >= 11 && abs(mcp->getPDG()) <= 16)  ) 
			{  
				// itself is a quark/lepton, jump over the intermediate paticle, search for the same quark/lepton, then ... 
				for ( unsigned int i=0; i<mcpa->getParents().size(); i++){
					MCParticle *mcpb = mcpa->getParents()[i];
					if ( mcpb->getPDG() == mcp->getPDG() ) {
						//printf( "here: %5d vs. %5d (%5d)\n", mcp->getPDG(),mcpb->getPDG(), mcpb->getParents()[0]->getPDG() );
						ret = FindParton(mcpb);
					}
				}
				if ( ret == NULL  && abs(mcp->getPDG()) <= 6 ) {
					for ( unsigned int i=0; i<mcpa->getParents().size(); i++ ) {
						MCParticle *mcpb = mcpa->getParents()[i];
						printf( "%5d vs. %5d, ", mcp->getPDG(),mcpb->getPDG() );
					}
					printf("\n");	
				}
			} 
			else 
			{  
				// itself is a hadron 
				// first; obtain partons from all grandparents
				map<MCParticle *, TVector3 > partonmap; partonmap.clear();
				map<MCParticle *, TVector3>::iterator it;
				for ( unsigned int i=0; i<mcpa->getParents().size(); i++){
					MCParticle *mcpb = mcpa->getParents()[i];
					MCParticle *mcpo = FindParton(mcpb);
					TVector3 v(mcpb->getMomentum());
					if ( mcpo != NULL) {
						it = partonmap.find(mcpo);
						if( it != partonmap.end() ) {
							v += it->second; // why add?
							//printf("the mother of %6d is %6d\n", (int)mcpb->getPDG(), (int)mcpo->getPDG());
						}
						partonmap[mcpo] = v;
					}
				}
				// then take the nearest one as the favorite parton
				double nearestangle = 2 * 3.14159;
				int selectedparton = -1, pdgid=0;
				MCParticle *pSelected = NULL;
				int j;
				TVector3 v2(mcp->getMomentum());
				for(j=0,it = partonmap.begin();it != partonmap.end(); it++,j++){
					if ( it->first ){
						double angle = v2.Angle(it->second);
						if(nearestangle > angle)
						{ 
							nearestangle = angle; selectedparton = j; pSelected = it->first; pdgid=(it->first)->getPDG();
						}
					}
				}
				ret = pSelected;
				if ( ret == NULL ) {
					//cout << partonmap.size() << " partons obtained." << endl;
					//cout << "Parton " << selectedparton << " selected. "<< nearestangle  << endl;
				}
				FreeAll(partonmap);
			}
		} 
		else
		{	
			ret = FindParton(mcpa);	
			if( 0 && ret==NULL){
				for (unsigned int k=0 ; k < mcpa->getParents().size() ; k++ ) {
					MCParticle 	* part    = mcpa->getParents()[k];
					printf("parent of %10d is %10d (%3d)\n", mcpa->getPDG(), part->getPDG(), (int)part->getParents().size() );
					for (unsigned int l=0 ; l < part->getParents().size() ; l++ ) {
						printf("parent %2d is %10d\n", l, part->getParents()[l]->getPDG() );
					}
				}
			}

		}	
	}
	else // More than one parent  
	{ 
		if ( abs(mcp->getPDG())<=6  ){
			MCParticle *mcpo2 = NULL;
			for( unsigned int i=0; i<mcp->getParents().size(); i++ ){
				MCParticle *mcpa = mcp->getParents()[i];
				if( mcp->getParents().size()==2 && abs(mcpa->getPDG()) == 11 ){
					mcpo2 = mcp;
					break;
				}else if( mcp->getPDG() == mcp->getPDG()  ){
					MCParticle *mcpo = FindParton(mcpa);
					mcpo2 = mcpo;
					break;
				}
			}
			ret=mcpo2;
		}
		else
		{	
			MCParticle *mcpo2 = NULL;
			for( unsigned int i=0; i<mcp->getParents().size(); i++ ){
				MCParticle *mcpa = mcp->getParents()[i];
				MCParticle *mcpo = FindParton(mcpa);
				if(mcpo2 && mcpo != mcpo2)
				{}
				//{cout << mcp->getPDG()<<" has Multiple parents other than intermidiates with different parton origin!" << endl;}
				mcpo2 = mcpo;
			}
			//
			ret = mcpo2;
		}
	}
	//
	return ret;
}
void LGFastJetClustering::analyze(const vector<fastjet::PseudoJet> & ee_kt_jets ) {

	/////// EnergyCorrelator /////////////////////////////

	// Initial clustering with ee_kt algorithm
	// fastjet::JetAlgorithm algorithm = fastjet::ee_kt_algorithm; 
	// fastjet::JetDefinition jetDef = fastjet::JetDefinition(algorithm);
	// fastjet::ClusterSequence clust_seq(input_particles,jetDef);
	// vector<fastjet::PseudoJet> ee_kt_jets  = fastjet::sorted_by_E(clust_seq.exclusive_jets(njet));

	int i =0; 
	for (unsigned int j = 0; j < ee_kt_jets.size(); j++) { // Two hardest jets per event
		if (ee_kt_jets[j].perp() > 0) {

			fastjet::PseudoJet myJet = ee_kt_jets[j];

			// various values of beta
			vector<double> betalist;
			betalist.push_back(0.1);
			betalist.push_back(0.2);
			betalist.push_back(0.5);
			betalist.push_back(1.0);
			betalist.push_back(1.5);
			betalist.push_back(2.0);

			// checking the two energy/angle modes
			vector<EnergyCorrelator::Measure> measurelist;
			//measurelist.push_back(EnergyCorrelator::pt_R);
			measurelist.push_back(EnergyCorrelator::E_theta); 

			vector<string> modename;
			//modename.push_back("pt_R");
			modename.push_back("E_theta");

			for (unsigned int M = 0; M < measurelist.size(); M++) {

				for (unsigned int B = 0; B < betalist.size(); B++) {
					double beta = betalist[B];

					EnergyCorrelator ECF0(0,beta,measurelist[M]);
					EnergyCorrelator ECF1(1,beta,measurelist[M]);
					EnergyCorrelator ECF2(2,beta,measurelist[M]);
					EnergyCorrelator ECF3(3,beta,measurelist[M]);
					EnergyCorrelator ECF4(4,beta,measurelist[M]);
					EnergyCorrelator ECF5(5,beta,measurelist[M]);

					_Ci[i++]=ECF1(myJet);
					_Ci[i++]=ECF2(myJet);
					_Ci[i++]=ECF3(myJet);
					_Ci[i++]=ECF4(myJet);
					_Ci[i++]=ECF5(myJet);
					printf("%7.3f %14.2f %14.2f %14.2f %14.2f %15.2f \n",beta,ECF1(myJet),ECF2(myJet),ECF3(myJet),ECF4(myJet),ECF5(myJet));
				}

				for (unsigned int B = 0; B < betalist.size(); B++) {
					double beta = betalist[B];

					EnergyCorrelatorRatio r0(0,beta,measurelist[M]);
					EnergyCorrelatorRatio r1(1,beta,measurelist[M]);
					EnergyCorrelatorRatio r2(2,beta,measurelist[M]);
					EnergyCorrelatorRatio r3(3,beta,measurelist[M]);
					EnergyCorrelatorRatio r4(4,beta,measurelist[M]);

					_Ci[i++]=r0(myJet);
					_Ci[i++]=r1(myJet);
					_Ci[i++]=r2(myJet);
					_Ci[i++]=r3(myJet);
					_Ci[i++]=r4(myJet);

					printf("%7.3f %14.4f %14.4f %14.4f %14.4f %15.4f \n",beta,r0(myJet),r1(myJet),r2(myJet),r3(myJet),r4(myJet));
				}

				for (unsigned int B = 0; B < betalist.size(); B++) {
					double beta = betalist[B];

					EnergyCorrelatorDoubleRatio C1(1,beta,measurelist[M]);
					EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M]);
					EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M]);
					EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M]);

					_Ci[i++]=C1(myJet);
					_Ci[i++]=C2(myJet);
					_Ci[i++]=C3(myJet);
					_Ci[i++]=C4(myJet);

					printf("%7.3f %14.6f %14.6f %14.6f %14.6f \n",beta,C1(myJet),C2(myJet),C3(myJet),C4(myJet));
				}
				cout << "-------------------------------------------------------------------------------------" << endl << endl;


			}
		}
	}
}



