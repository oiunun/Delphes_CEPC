#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <UTIL/PIDHandler.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "CheckJetClusteringProcessor.h"


using namespace std;

CheckJetClusteringProcessor aCheckJetClusteringProcessor ;


CheckJetClusteringProcessor::CheckJetClusteringProcessor() : Processor("CheckJetClusteringProcessor") {

	// modify processor description
	_description = "CheckJetClusteringProcessor does whatever it does ..." ;

	// register steering parameters: name, description, class-variable, default value
	registerInputCollection( LCIO::MCPARTICLE,
			"InputMCParticlesCollection" , 
			"Name of the MCParticle collection",
			_colMCP ,
			std::string("MCParticle") ) ;


	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputArborPFOsCollection" , 
			"Name of the PFOs collection"  ,
			_colPFOs ,
			std::string("ArborPFOs") ) ;


	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputPandoraPFOsCollection" , 
			"Name of the PFOs collection"  ,
			_colPFOs ,
			std::string("PandoraPFOs") ) ;

	registerInputCollection( LCIO::LCRELATION,
			"InputMCTruthLinkCollection" , 
			"Name of the MCTruthLink collection"  ,
			_colMCTLJet ,
			std::string("RecoMCTruthLink") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputJetsCollection",
			"Name of the jets collection",
			_colJets,
			std::string("RefinedJets") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputPDRJetsCollection",
			"Name of the jets collection",
			_panJets,
			std::string("panJets") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputJetsCollectionFaJ",
			"Name of the jets collection",
			_colJetF,
			std::string("RefinedJets") );

	registerInputCollection( LCIO::LCRELATION,
			"InputMCTruthLinkCollectionGen" , 
			"Name of the MCTruthLink collection Gen"  ,
			_colMCTLGen ,
			std::string("RecoMCTruthMap") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputJetsCollectionGen",
			"Name of the jets collection",
			_colGenJ,
			std::string("GenJets") );


	registerProcessorParameter("MatchMC",          "1: match MC; 0: not match",  m_matchmc,         0);

	registerProcessorParameter("SaveHist",         "1: save histogram; 0: not",  m_savehis,         0);

	registerProcessorParameter("DEBUG",            "Verbose level",              m_debug,           0);

	registerProcessorParameter("ShowMC",           "print mc  or not",           m_showmc,          0);
	registerProcessorParameter("FastJet",          "pure fastjet ",              m_fastjet,         0);

	registerProcessorParameter("ECM",              "Energy of C.M. ",            m_ECM,         250.0);

	registerProcessorParameter("EnergyCut",        "Threshold of all particles", m_EnergyCut,           0.5);

	registerProcessorParameter("LepEnergy",        "Threshold of all particles", m_LepEnergyThreshold,  0.5);

}



void CheckJetClusteringProcessor::init() { 

	streamlog_out(DEBUG) << "   init called  " << std::endl ;

	// usually a good idea to
	printParameters() ;
	m_ecms             = TLorentzVector(0,0,0,m_ECM);

	//
	m_FSNames[0]="EXC4_0";
	for (unsigned int i = 0; i < 1; i++){
		m_CutPass [i]=0;
		m_Checking[i]=0;
		for (int j = 0; j < 20; j++) m_cutpass[i][j]=0;
		if (m_FSNames[i] != ""){
			// a tree for reconstructed information
			char* ntFullName = new char[100];
			sprintf(ntFullName,"nt%s","CheckJet");
			TTree* outputTree = new TTree(ntFullName,"");
			outputTree->SetAutoSave(32*1024*1024); 
			NTupleHelper* nt = new NTupleHelper( outputTree, m_ecms ); 
			delete ntFullName;
			// a tree for generated information
			NTupleHelper* ntgen = NULL;
			string ntGenName("ntGEN");  ntGenName += m_FSNames[i];
			// create an FSFinfo object
			FSInfo* fsInfo = new FSInfo(m_FSNames[i],nt,ntgen);
			m_FSInfoVector.push_back(fsInfo);

		}
	}

	_nRun = 0 ;
	_nEvt = 0 ;
}

void CheckJetClusteringProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void CheckJetClusteringProcessor::processEvent( LCEvent * evt ) { 
	_nEvt ++ ;
	_FD   =-1.;
	_Pmax = 0 ;
	_Emax = 0 ;
	if( m_debug>1 ) printf("\n");
	if( m_debug>1 || _nEvt%100==0 ) printf("Run & event =  %8d %8d\n",_nRun, _nEvt);
	m_HiggsBoostVector = TVector3(0,0,0); 
	//
	LCCollection    *colMC          = evt->getCollection(_colMCP );
	LCCollection    *colPFOs        = NULL; 
	LCCollection    *colPDRs        = NULL; 
	LCCollection    *colJets        = NULL; 
	LCCollection    *panJets        = NULL; 
	LCCollection    *colJetF        = NULL; 
	LCCollection    *colMCTLJet     = NULL;
	LCCollection    *colGenJ        = NULL; 
	LCCollection    *colMCTLGen     = NULL;
	//	
	if( _colPFOs != "" ) colPFOs = evt->getCollection( _colPFOs );
	if( _colPDRs != "" ) colPDRs = evt->getCollection( _colPDRs );
	//	
	//	
	if( _colJets != ""){
		colJets  = evt->getCollection(_colJets);
		panJets  = evt->getCollection(_panJets);
		colGenJ  = evt->getCollection(_colGenJ);
		colJetF  = evt->getCollection(_colJetF);
		if( m_debug>1 ) { 
			printf("No. of elements in  JetCol = %3d\n", colJets->getNumberOfElements()    );
			printf("No. of elements in  JetPDR = %3d\n", panJets->getNumberOfElements()    );
			printf("No. of elements in  GenJet = %3d\n", colGenJ->getNumberOfElements()    );
			printf("No. of elements in  FajJet = %3d\n", colGenJ->getNumberOfElements()    );
		}
	}
	//	
	if( m_matchmc>0 ) {
		colMCTLJet = evt->getCollection( _colMCTLJet );
		colMCTLGen = evt->getCollection( _colMCTLGen );
	}
	//
	m_mcTruthHelper = new MCTruthHelper(colMC);
	if( m_showmc ) m_mcTruthHelper->printInformation(m_showmc);

	_nb = m_mcTruthHelper->getnBquark();
	_nc = m_mcTruthHelper->getnCquark();
	_nl = m_mcTruthHelper->getnLquark();
	//	
	_nhgg = m_mcTruthHelper->nHiggsDaughters(21);
	_nhpp = m_mcTruthHelper->nHiggsDaughters(22);
	_nhzz = m_mcTruthHelper->nHiggsDaughters(23);
	_nhww = m_mcTruthHelper->nHiggsDaughters(24);
	_nhee = m_mcTruthHelper->nHiggsDaughters(11);
	_nhmm = m_mcTruthHelper->nHiggsDaughters(13);
	_nhtt = m_mcTruthHelper->nHiggsDaughters(15);
	_nhuu = m_mcTruthHelper->nHiggsDaughters( 1);
	_nhdd = m_mcTruthHelper->nHiggsDaughters( 2); 
	_nhss = m_mcTruthHelper->nHiggsDaughters( 3);
	_nhcc = m_mcTruthHelper->nHiggsDaughters( 4);
	_nhbb = m_mcTruthHelper->nHiggsDaughters( 5);
	// _nhfs represent the  higgs decay final state
	// 1-6: dd, uu, ss, cc, bb, tt
	// 7: ee; 8: mumu; 9 tautau; 10: gluongluon;
	// 11: gammagamma, 12: gamma Z; 13: ZZ, 14: WW
	_nhfs = m_mcTruthHelper->nHiggsFinalState( );
	//
	BuildParticleList( colMC, colPFOs, colPDRs, colJets, colJetF, colMCTLJet, colGenJ, colMCTLGen, panJets );
	//
	//********************************************************************
	//
	//   LOOP OVER FINAL STATES
	//
	//********************************************************************
	//
	for (unsigned int ifs = 0; ifs < m_FSInfoVector.size(); ifs++){
      if( JetList.size()<4 || GenList.size()<4 || PanList.size()<4) continue;
		m_cutpass[ifs][0]++;

		// get information about this final state
		FSInfo* fsinfo                = m_FSInfoVector[ifs];
		//fsinfo->Print();
		NTupleHelper* NT              = fsinfo->NT();
		vector<string> particleNames  = fsinfo->particleNames();
		vector<int>    particleStatus = fsinfo->particleStatus();
		m_cutpass[ifs][1]++;

		//********************************************************************
		//
		//   PUT TOGETHER ALL COMBINATIONS OF PARTICLES FOR THIS FINAL STATE
		//
		//********************************************************************
		vector< vector< FSParticle* > > pCombos; pCombos.clear();
		vector< FSParticle* > combo1; combo1.clear();
		vector< FSParticle* > combo2; combo2.clear();
		vector< FSParticle* > combo3; combo3.clear();
		vector< FSParticle* > combo4; combo4.clear();
		for ( unsigned int ipl = 0; ipl < JetList.size(); ipl++ ){
			combo1.push_back(JetList[ipl]); 
			combo2.push_back(GenList[ipl]); 
			combo3.push_back(FaJList[ipl]); 
			combo4.push_back(PanList[ipl]); 
		}
		pCombos.push_back(combo1);
		pCombos.push_back(combo2);
		pCombos.push_back(combo3);
		pCombos.push_back(combo4);
		m_cutpass[ifs][2]++;
		//	
		if( pCombos.size() < 3 ) continue;
		m_cutpass[ifs][3]++;
		NT->fillDouble("Run",         _nRun);
		NT->fillDouble("Event",       _nEvt);
		NT->fillDouble("ntrks",       ntrks);
		NT->fillDouble("nclus",       nclus);
		NT->fillDouble("nPFOs",       nPFOs);
		NT->fillDouble("nPDRs",       nPDRs);
		NT->fillDouble("Pmax",        _Pmax);
		NT->fillDouble("Emax",        _Emax);
		NT->fillDouble("njets",       numberJets);
		NT->fillDouble("VisEn",       pVis.E());
		NT->fillDouble("VisPx",       pVis.X());
		NT->fillDouble("VisPy",       pVis.Y());
		NT->fillDouble("VisPz",       pVis.Z());
		NT->fillDouble("VisMass",     pVis.M());
		NT->fillDouble("Fuly34",      _y34_jet);
		NT->fillDouble("Fuly45",      _y45_jet);
		NT->fillDouble("Geny34",      _y34_gen);
		NT->fillDouble("Geny45",      _y45_gen);
		if(m_savemc)NT->fillMCTruth(m_mcTruthHelper);
		NT->fillDouble("nhfs" ,    _nhfs  );
		NT->fillDouble("VisEnMC",  m_mcTruthHelper->MCTotalEnergy()   );
		NT->fillDouble("MisEnMC",  m_mcTruthHelper->MCMissingEnergy() );
		NT->fillDouble("HiggsMass",m_mcTruthHelper->getHiggsP4(colMC).M() );
		NT->fillDouble("ZMass",    m_mcTruthHelper->getZP4(colMC).M() );
		NT->fillDouble("E_nu_H",    m_mcTruthHelper->getE_nu_H()  );
		NT->fillDouble("E_nu_Z",    m_mcTruthHelper->getE_nu_Z()  );
		NT->fill4Momentum("H", m_mcTruthHelper->getHiggsP4(colMC));
		NT->fill4Momentum("Z", m_mcTruthHelper->getZP4(colMC));
		//printf("H = %10.2f, Z =%10.2f\n",  m_mcTruthHelper->getHiggsP4(colMC).M(), m_mcTruthHelper->getZP4(colMC).M() ); 
		//
		m_cutpass[ifs][4]++;
		//********************************************************************
		//
		//   CUT ON THE TOTAL ENERGY AND MOMENTUM TO SAVE TIME
		//
		//********************************************************************
		vector<FSParticle*> combo = pCombos[0];

		totalP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
		for (unsigned int t = 0; t < combo.size(); t++){
			if(!combo[t]->missed() ) totalP4 += combo[t]->rawFourMomentum();
		}
		double missingMass2  = (totalP4-m_ecms).M2();
		NT->fillDouble("FulMissingMass2",   missingMass2 );
		NT->fillDouble("FulTotalP",         totalP4.Rho());
		NT->fillDouble("FulTotalEnergy",    totalP4.E()  );
		NT->fillDouble("FulTotalPx",        totalP4.Px() );
		NT->fillDouble("FulTotalPy",        totalP4.Py() );
		NT->fillDouble("FulTotalPz",        totalP4.Pz() );

		m_cutpass[ifs][5]++;
		m_cutpass[ifs][6]++;
		m_cutpass[ifs][7]++;
		m_cutpass[ifs][8]++;
		m_cutpass[ifs][9]++;
		//********************************************************************
		//
		//   RECORD INFORMATION
		//
		//********************************************************************
		vector<TLorentzVector> rawp4list;
		// full jet
		rawp4list.clear();  
		char index[30];
		for (unsigned int j = 0; j < combo.size(); j++){
			FSParticle* fsp = combo[j];
			rawp4list.push_back(fsp->rawFourMomentum());
			NT->fillJet(fsp, j+1,"FulJet", 1 , m_kappa);
		}
		//
		NT->fill4Momentum("idx_raw","Ful_", rawp4list, rawp4list.size());
		if( rawp4list.size()>=2 ){
			for(unsigned int ki=0; ki<rawp4list.size()-1; ki++){
				for(unsigned int kj=ki+1; kj<rawp4list.size(); kj++){
					sprintf(index ,"FulMass%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,(rawp4list[ki]+rawp4list[kj]).M());
					sprintf(index ,"FulCosTheta%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,cos(rawp4list[ki].Angle(rawp4list[kj].Vect())));
				}
			}
		}
		// genjet 
		combo = pCombos[1];
		rawp4list.clear();  
		totalP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
		for (unsigned int t = 0; t < combo.size(); t++){
			if(!combo[t]->missed() ) totalP4 += combo[t]->rawFourMomentum();
		}
		missingMass2  = (totalP4-m_ecms).M2();
		NT->fillDouble("GenMissingMass2",   missingMass2 );
		NT->fillDouble("GenTotalP",         totalP4.Rho());
		NT->fillDouble("GenTotalEnergy",    totalP4.E()  );
		NT->fillDouble("GenTotalPx",        totalP4.Px() );
		NT->fillDouble("GenTotalPy",        totalP4.Py() );
		NT->fillDouble("GenTotalPz",        totalP4.Pz() );
		rawp4list.clear();  
		for (unsigned int j = 0; j < combo.size(); j++){
			FSParticle* fsp = combo[j];
			rawp4list.push_back(fsp->rawFourMomentum());
			NT->fillJet(fsp, j+1,"GenJet", 0 , m_kappa);
		}
		//
		NT->fill4Momentum("idx_raw","Gen_", rawp4list, rawp4list.size());
		if( rawp4list.size()>=2 ){
			for(unsigned int ki=0; ki<rawp4list.size()-1; ki++){
				for(unsigned int kj=ki+1; kj<rawp4list.size(); kj++){
					sprintf(index ,"GenMass%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,(rawp4list[ki]+rawp4list[kj]).M());
					sprintf(index ,"GenCosTheta%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,cos(rawp4list[ki].Angle(rawp4list[kj].Vect())));
				}
			}
		}
		//fast jet 
		combo = pCombos[2];
		totalP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
		for (unsigned int t = 0; t < combo.size(); t++){
			if(!combo[t]->missed() ) totalP4 += combo[t]->rawFourMomentum();
		}
		missingMass2  = (totalP4-m_ecms).M2();
		NT->fillDouble("FaJMissingMass2",   missingMass2 );
		NT->fillDouble("FaJTotalP",         totalP4.Rho());
		NT->fillDouble("FaJTotalEnergy",    totalP4.E()  );
		NT->fillDouble("FaJTotalPx",        totalP4.Px() );
		NT->fillDouble("FaJTotalPy",        totalP4.Py() );
		NT->fillDouble("FaJTotalPz",        totalP4.Pz() );
		rawp4list.clear();  
		for (unsigned int j = 0; j < combo.size(); j++){
			FSParticle* fsp = combo[j];
			rawp4list.push_back(fsp->rawFourMomentum());
			NT->fillJet(fsp, j+1,"FaJJet", 1 , m_kappa);
		}
		//
		NT->fill4Momentum("idx_raw","FaJ_", rawp4list, rawp4list.size());
		if( rawp4list.size()>=2 ){
			for(unsigned int ki=0; ki<rawp4list.size()-1; ki++){
				for(unsigned int kj=ki+1; kj<rawp4list.size(); kj++){
					sprintf(index ,"FaJMass%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,(rawp4list[ki]+rawp4list[kj]).M());
					sprintf(index ,"FaJCosTheta%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,cos(rawp4list[ki].Angle(rawp4list[kj].Vect())));
				}
			}
		}
		// pandora jet 
		combo = pCombos[3];
		totalP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
		for (unsigned int t = 0; t < combo.size(); t++){
			if(!combo[t]->missed() ) totalP4 += combo[t]->rawFourMomentum();
		}
		missingMass2  = (totalP4-m_ecms).M2();
		NT->fillDouble("PDRMissingMass2",   missingMass2 );
		NT->fillDouble("PDRTotalP",         totalP4.Rho());
		NT->fillDouble("PDRTotalEnergy",    totalP4.E()  );
		NT->fillDouble("PDRTotalPx",        totalP4.Px() );
		NT->fillDouble("PDRTotalPy",        totalP4.Py() );
		NT->fillDouble("PDRTotalPz",        totalP4.Pz() );
		rawp4list.clear();  
		for (unsigned int j = 0; j < combo.size(); j++){
			FSParticle* fsp = combo[j];
			rawp4list.push_back(fsp->rawFourMomentum());
			NT->fillJet(fsp, j+1,"PDRJet", 1 , m_kappa);
		}
		//
		NT->fill4Momentum("idx_raw","PDR_", rawp4list, rawp4list.size());
		if( rawp4list.size()>=2 ){
			for(unsigned int ki=0; ki<rawp4list.size()-1; ki++){
				for(unsigned int kj=ki+1; kj<rawp4list.size(); kj++){
					sprintf(index ,"FaJMass%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,(rawp4list[ki]+rawp4list[kj]).M());
					sprintf(index ,"FaJCosTheta%d%d",ki+1,kj+1);
					NT->fillDouble((string)index,cos(rawp4list[ki].Angle(rawp4list[kj].Vect())));
				}
			}
		}
		// write the tree
		m_CutPass[ifs]++;
		NT->write();
		//
		FreeAll(rawp4list);
	}
	//
	CleanEvt();
}

void CheckJetClusteringProcessor::check( LCEvent * evt ) 
{ 

}


void CheckJetClusteringProcessor::end()
{
	char message[10][40]={
		"Input Number of   Evts ",//  0 
		"nChrg  protection  Cut ",//  1
		"nCombo uplimit protect ",//  2
		"nCombo greater than  0 ",//  3
		"Before  E and P    Cut ",//  4
		"Missing Energy     Cut ",//  5
		"Missing Momentum   Cut ",//  6
		"Missing Mass       Cut ",//  7
		"Raw   4-Momentum   Cut ",//  8
		"No  of Filling Entries " //  9
	};	

	printf("\n");
	cout << "Total events = " << _nEvt << endl;
	for(unsigned int i=0; i<m_FSInfoVector.size(); i++){
		printf("\n");
		printf("%27s\n", (m_FSInfoVector[i]->FSName()).data());
		for(int j=0; j<10; j++){
			if     (j==0)printf("%3d) %30s: %10d %10.2f%%\n",j,message[j],(int)m_cutpass[0][j],m_cutpass[0][j]/m_cutpass[0][j-0]*100);
			else if(j< 2)printf("%3d) %30s: %10d %10.2f%%\n",j,message[j],(int)m_cutpass[0][j],m_cutpass[0][j]/m_cutpass[0][j-1]*100);
			else if(j==2)printf("%3d) %30s: %10d %10.2f%%\n",j,message[j],(int)m_cutpass[i][j],m_cutpass[i][j]/m_cutpass[0][j-1]*100);
			else         printf("%3d) %30s: %10d %10.2f%%\n",j,message[j],(int)m_cutpass[i][j],m_cutpass[i][j]/m_cutpass[i][j-1]*100);
		}
	}

	printf("\n");
	for(unsigned int i=0; i<m_FSInfoVector.size(); i++){
		printf("%20s fill number :%10d\n", (m_FSInfoVector[i]->FSName()).data(), m_CutPass[i]);
	}

	printf("\n\n\n");

	FreeDelAll(m_FSInfoVector);	
}

void
CheckJetClusteringProcessor::BuildParticleList(
		LCCollection    *colMC  , LCCollection    *colPFO , LCCollection    *colPDR , 
		LCCollection   *colJets , LCCollection    *colJetF, LCCollection    *colMCTLJet,
		LCCollection   *colGenJ , LCCollection    *colMCTLGen, LCCollection    *panJets
		){

	if( colPFO ){
		_Pmax = 0 ;
		_Emax = 0 ;
		nPFOs = colPFO->getNumberOfElements();
		raw_PFOs.clear(); 
		//
		for (int i=0;i<nPFOs;i++) {
			ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
			Double_t energy = recPart->getEnergy();
			if ( energy<m_EnergyCut ) continue; 
			raw_PFOs.push_back(recPart); 
			TrackVec   tckvec = recPart->getTracks();
			ClusterVec cluvec = recPart->getClusters();
			int ntrk = tckvec.size();
			int nclu = cluvec.size();
			//
			Double_t charge = recPart->getCharge();
			TVector3 momentum = TVector3(recPart->getMomentum());
			Double_t momentumMagnitude = momentum.Mag();
			ntrks += ntrk;
			if(fabs(charge)>0.1)nclus += nclu;
			if( energy   > _Emax ) _Emax = energy  ;
			if( momentumMagnitude  > _Pmax ) _Pmax = momentumMagnitude;
			TLorentzVector p4 = TLorentzVector(momentum,energy);
			if ( energy>m_EnergyCut   && energy<m_ECM*0.6) {
				pVis += p4;
			}
		}
		nPFOs = raw_PFOs.size(); 
	}
	if( colPDR ){
		nPDRs = colPDR->getNumberOfElements();
		raw_PDRs.clear(); 
		//
		for (int i=0;i<nPDRs;i++) {
			ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPDR->getElementAt(i));
			Double_t energy = recPart->getEnergy();
			if ( energy<m_EnergyCut ) continue; 
			raw_PDRs.push_back(recPart); 
		}
		nPDRs = raw_PDRs.size(); 
	}


	vector<MCParticleImpl*>  partonList = m_mcTruthHelper->partonList();
	_y34_jet =  -1.0;
	_y45_jet =  -1.0;
	_y34_gen =  -1.0;
	_y45_gen =  -1.0;
	_y34_faj =  -1.0;
	_y45_faj =  -1.0;
	if( colJets ){
		numberJets = colJets ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of jets in Col= %3d\n",numberJets ) ;
		}
		if ( numberJets > 0) {
			double btag=-1, ctag=-1, bctag=-1, cat=-1, flavor=0;
			for (int j=0 ; j < numberJets ; j++ ) {
				ReconstructedParticle *pJet = dynamic_cast<ReconstructedParticle *>( colJets->getElementAt (j) );
				TVector3 momentum = TVector3(pJet->getMomentum());
				JetList .push_back(new FSParticle(pJet, colMCTLJet, partonList,  "jet", btag, ctag, bctag, cat, flavor));
			}
			if(m_fastjet==0){
				PIDHandler pidh(colJets);
				ReconstructedParticle *Jet0 = dynamic_cast<ReconstructedParticle *>(colJets->getElementAt(0));
				Int_t algo_y = pidh.getAlgorithmID("yth");
				const ParticleID & ythID = pidh.getParticleID(Jet0, algo_y);
				FloatVec params_y = ythID.getParameters();
				_y34_jet = params_y[pidh.getParameterIndex(algo_y, "y34")];
				_y45_jet = params_y[pidh.getParameterIndex(algo_y, "y45")];
				_y34_jet =  _y34_jet>1.0001 ?  _y34_jet-100.0 : _y34_jet;  
				_y45_jet =  _y45_jet>1.0001 ?  _y45_jet-100.0 : _y45_jet;
			}else{	
				_y34_jet =  colJets->parameters().getFloatVal( "y34" );
				_y45_jet =  colJets->parameters().getFloatVal( "y45" );
			}
		}
	}

	if( panJets ){
		numberJets = panJets ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of jets in PDR= %3d\n",numberJets ) ;
		}
		if ( numberJets > 0) {
			double btag=-1, ctag=-1, bctag=-1, cat=-1, flavor=0;
			for (int j=0 ; j < numberJets ; j++ ) {
				ReconstructedParticle *pJet = dynamic_cast<ReconstructedParticle *>( panJets->getElementAt (j) );
				TVector3 momentum = TVector3(pJet->getMomentum());
				PanList .push_back(new FSParticle(pJet, colMCTLGen, partonList,  "jet", btag, ctag, bctag, cat, flavor));
			}
			_y34_pan =  panJets->parameters().getFloatVal( "y34" );
			_y45_pan =  panJets->parameters().getFloatVal( "y45" );
		}

	}


	if( colGenJ ){
		numberJets = colGenJ ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of jets in Col= %3d\n",numberJets ) ;
		}
		if ( numberJets > 0) {
			double btag=-1, ctag=-1, bctag=-1, cat=-1, flavor=0;
			for (int j=0 ; j < numberJets ; j++ ) {
				ReconstructedParticle *pJet = dynamic_cast<ReconstructedParticle *>( colGenJ->getElementAt (j) );
				TVector3 momentum = TVector3(pJet->getMomentum());
				GenList .push_back(new FSParticle(pJet, colMCTLGen, partonList,  "jet", btag, ctag, bctag, cat, flavor));
			}
			_y34_gen =  colGenJ->parameters().getFloatVal( "y34" );
			_y45_gen =  colGenJ->parameters().getFloatVal( "y45" );
		}
	}

	if( colJetF ){
		numberJets = colJetF ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of jets in Col= %3d\n",numberJets ) ;
		}
		if ( numberJets > 0) {
			double btag=-1, ctag=-1, bctag=-1, cat=-1, flavor=0;
			for (int j=0 ; j < numberJets ; j++ ) {
				ReconstructedParticle *pJet = dynamic_cast<ReconstructedParticle *>( colJetF->getElementAt (j) );
				TVector3 momentum = TVector3(pJet->getMomentum());
				FaJList .push_back(new FSParticle(pJet, colMCTLJet, partonList,  "jet", btag, ctag, bctag, cat, flavor));
			}
			_y34_faj  =  colJetF->parameters().getFloatVal( "y34" );
			_y45_faj  =  colJetF->parameters().getFloatVal( "y45" );
		}
	}

}	

void
CheckJetClusteringProcessor::CleanEvt()
{
	FreeDelAll(  JetList      );
	FreeDelAll(  PanList      );
	FreeDelAll(  GenList      );
	FreeDelAll(  FaJList      );
	FreeDelAll( ParticleTrash );
	FreeAll   ( DictPList     );
	delete m_mcTruthHelper; m_mcTruthHelper=NULL;
}

bool
CheckJetClusteringProcessor::checkCombination(const vector<FSParticle*>& combo, bool complete, bool inclusive){

	if (inclusive) return true;
	// if the combination isn't yet complete, just check to make sure there
	//  is no excess energy outside of the energy tolerance
	if (!complete){
		double totalEnergy = 0.0;
		for (unsigned int i = 0; i < combo.size(); i++){
			if(!combo[i]->missed()) totalEnergy += combo[i]->rawFourMomentum().E();
		}

		//double excessEnergy = totalEnergy - m_ecms.E();
		//if (excessEnergy > m_energyTolerance) return false;
		return true;

	}

	// if the combination is complete, calculate the total energy and momentum
	TLorentzVector pTotal(0.0, 0.0, 0.0, 0.0);
	for (unsigned int i = 0; i < combo.size(); i++){
		if(!combo[i]->missed())	pTotal += combo[i]->rawFourMomentum();
	}

	//HepLorentzVector pMissing = m_ecms - pTotal;

	// if the combination is complete and exclusive, check the energy
	//   and momentum balance
	//if (fabs(pMissing.e())          > m_energyTolerance  ) return false;
	//if (fabs(pMissing.vect().mag()) > m_momentumTolerance) return false;
	return true;
}
