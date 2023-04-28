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
#include "FSClasserProcessor.h"


using namespace std;

FSClasserProcessor aFSClasserProcessor ;


FSClasserProcessor::FSClasserProcessor() : Processor("FSClasserProcessor") {

	// modify processor description
	_description = "FSClasserProcessor does whatever it does ..." ;

	// register steering parameters: name, description, class-variable, default value
	registerInputCollection( LCIO::MCPARTICLE,
			"InputMCParticlesCollection" , 
			"Name of the MCParticle collection",
			_colMCP ,
			std::string("MCParticle") ) ;

	registerInputCollection( LCIO::LCRELATION,
			"InputMCTruthLinkCollection" , 
			"Name of the MCTruthLink collection"  ,
			_colMCTL ,
			std::string("RecoMCTruthLink") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputPandoraPFOsCollection" , 
			"Name of the PFOs collection"  ,
			_colPFOs ,
			std::string("ArborPFOs") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputIsoLepsCollection" , 
			"Name of the Isolated Leptons collection"  ,
			_colIsoLeps ,
			std::string("") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"InputIsoZLepsCollection" , 
			"Name of the Isolated Leptons collection from Z decay"  ,
			_colIsoZLeps ,
			std::string("") ) ;

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputTauCollection",
			"Name of the Tau collection",
			_colTaus,
			std::string("") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputJetsCollection",
			"Name of the 2 jets collection",
			_col2Jets,
			std::string("") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputJets4Collection",
			"Name of the 4 jets collection",
			_col4Jets,
			std::string("") );


	registerProcessorParameter("MatchMC",          "1: match MC; 0: not match",  m_matchmc,         0);

	registerProcessorParameter("SaveHist",         "1: save histogram; 0: not",  m_savehis,         0);

	registerProcessorParameter("SaveMC",           "1: save MC; 0: not save",    m_savemc,          0);

	registerProcessorParameter("MakePlots",        "Make some plots for check",  m_makeplots,       0);

	registerProcessorParameter("DEBUG",            "Verbose level",              m_debug,           0);

	registerProcessorParameter("kmfit",            "kinematic fit or not",       m_kmfit,           0);

	registerProcessorParameter("useCalo",          "use calo info for PID",      m_useCalo,         0);

	registerProcessorParameter("Luxury",           "save extra tracks/photons",  m_luxury,          0);

	registerProcessorParameter("ShowMC",           "print mc  or not",           m_showmc,          0);

	registerProcessorParameter("FastOrFull",       "Fass or Full simul.",        m_full,            1);

	registerProcessorParameter("ECM",              "Energy of C.M. ",            m_ECM,         250.0);

	registerProcessorParameter("Weight",           "Weight of evts, d=1",        m_wate,          1.0);

	registerProcessorParameter("Kappa",            "jet charge factor",          m_kappa,         0.3);

	registerProcessorParameter("TagFlavor",        "use flavor Tag??",           m_TagFlavor,       0);

	registerProcessorParameter("EventShape",       "Cal. evt-shape parameters",  m_EventShape,      0);
	
	registerProcessorParameter("LinearSphericity", "Cal. linear/sphericity",     m_LinearSphericity,0);

   //the final states of recoil side to calculate the boost vector, such as uu of Z(uu)H. NULL for no boost taken into account	
	registerProcessorParameter("RecFS4Boost",      "final state of Higgs recoil side", m_boostFS,string(""));
	
	registerProcessorParameter("EnergyCut",        "Threshold of all particles", m_EnergyCut,           0.5);

	registerProcessorParameter("LepEnergy",        "Threshold of all particles", m_LepEnergyThreshold,  0.5);

	for(unsigned int i=0; i<MAXFS; i++){
		char fsname[100];
		sprintf(fsname,"FS%3.3d",i);
		registerProcessorParameter(fsname, "Description of final state", m_FSNames[i], std::string(""));
	}

	for(unsigned int i=0; i<MAXFS; i++){
		char fsname[100];
		sprintf(fsname,"FSCut%3.3d",i);
		registerProcessorParameter(fsname, "Cuts of the final state"   , m_FSCuts[i], vector<string>());
	}
}



void FSClasserProcessor::init() { 
  
	streamlog_out(DEBUG) << "   init called  " << std::endl ;

	// usually a good idea to
	printParameters() ;
	double   pz        = m_ECM*gRandom->Gaus(0., 0.0017);
	m_ecms             = TLorentzVector(0, 0, pz, m_ECM+pz); 
	m_findEp           = false;
	m_findEm           = false;
	m_findMup          = false;
	m_findMum          = false;
	m_findTaup         = false;
	m_findTaum         = false;
	m_findJets         = false;
	m_findPhoton       = false;
	//
	for (unsigned int i = 0; i < MAXFS; i++){
		m_CutPass [i]=0;
		m_Checking[i]=0;
		for (int j = 0; j < 20; j++) m_cutpass[i][j]=0;
		if (m_FSNames[i] != ""){
			// a tree for reconstructed information
			char* ntFullName = new char[100];
			sprintf(ntFullName,"nt%s",m_FSNames[i].data());
			TTree* outputTree = new TTree(ntFullName,"");
			outputTree->SetAutoSave(32*1024*1024); 
			NTupleHelper* nt = new NTupleHelper( outputTree, m_ecms ); 
			delete ntFullName;
			// a tree for generated information
			string ntGenName("ntGEN");  ntGenName += m_FSNames[i];
			NTupleHelper* ntgen = NULL;
			if (ntGenName.find("EXC") != string::npos){
				ntGenName.erase(ntGenName.find("EXC"),3);
				string ntgenfullname("");  ntgenfullname += ntGenName;
				char* ntGenFullName = new char[100];
				sprintf(ntGenFullName,"%s",ntgenfullname.data());
				TTree* outputTree = new TTree(ntGenFullName, "");
				outputTree->SetAutoSave(32*1024*1024); 
				ntgen =  new NTupleHelper( outputTree, m_ecms ); 
				delete ntGenFullName;
			}
			// create an FSFinfo object
			FSInfo* fsInfo = new FSInfo(m_FSNames[i],nt,ntgen);
			m_FSInfoVector.push_back(fsInfo);

			// check to see which particles we need to reconstruct
			if ( fsInfo->hasParticle("e+"     ) ) m_findEp      = true;
			if ( fsInfo->hasParticle("e-"     ) ) m_findEm      = true;
			if ( fsInfo->hasParticle("mu+"    ) ) m_findMup     = true;
			if ( fsInfo->hasParticle("mu-"    ) ) m_findMum     = true;
			if ( fsInfo->hasParticle("tau+"   ) ) m_findTaup    = true;
			if ( fsInfo->hasParticle("tau-"   ) ) m_findTaum    = true;
			if ( fsInfo->hasParticle("gamma"  ) ) m_findPhoton  = true;
			if ( fsInfo->hasParticle( "jet"   ) ) m_findJets    = true;
		}
	}

	for (unsigned int i = 0; i < m_FSInfoVector.size(); i++){
		FSInfo* fs = m_FSInfoVector[i];
		int n = fs->nMissingParticles();
		if( n>1) {
			printf("\nChannel %3d: %s ",i, fs->FSName().data());
			printf("has %3d missing particles, please check ... \n",n);
			fs->Print();
			exit(1);        
		}               
	} 

	for (unsigned int i = 0; i < m_FSInfoVector.size(); i++){
		FSInfo* fs = m_FSInfoVector[i];
		vector<string> particleNames = fs->particleNames();
		printf("Channel %3d: %s\n",i, fs->FSName().data());
		fs->Print();                    
	} 

	for (unsigned int i = 0; i < MAXFS; i++){
		if ( m_FSCuts[i].size()!=5 ) continue;
		if ( m_FSCuts[i][0] != ""  ){
			FSCut* fscut = new FSCut(m_FSCuts[i]);
			printf("Cut of %15s, for comb=%s, type = %s, cut = %15.3f %15.3f\n", 
					fscut->FSName().data(), fscut->submodeName().data(), 
					fscut->cutType().data(), fscut->lowCut(), fscut->highCut() );
			bool found = false;
			for (unsigned int ifs = 0; ifs < m_FSInfoVector.size(); ifs++){
				if (m_FSInfoVector[ifs]->FSName() == fscut->FSName()){
					m_FSInfoVector[ifs]->addFSCut(fscut);
					found = true;
					break;
				}
			}
			if (!found){
				cout << "FSClasser ERROR: could not find a final state for the " << endl;
				cout << "FSCut with these arguments: " << endl;
				cout << m_FSCuts[i][0] << endl;
				cout << "  and final state = " << fscut->FSName() << endl;
				exit(0);
			}
		}

	}
	_nRun = 0 ;
	_nEvt = 0 ;
	for(int i=0; i<20; i++) hValid[i] = NULL;

	//if( m_full==0 ) { m_savehis=0; m_makeplots=0; }
	if( m_savehis ){
		//
		char* ntFullName = new char[100];
		sprintf(ntFullName,"ntValidation");
		TTree* outputTree = new TTree(ntFullName,"");
		outputTree->SetAutoSave(32*1024*1024); 
		m_ntVal = new NTupleHelper( outputTree, m_ecms ); 
		delete ntFullName;
		//
		hValid[ 0] = new TH1D("Valid_00", "Total Energy",100,   0.0, 300.0);
		hValid[ 1] = new TH1D("Valid_01", "Total Px    ",100,-100.0, 100.0);
		hValid[ 2] = new TH1D("Valid_02", "Total Py    ",100,-100.0, 100.0);
		hValid[ 3] = new TH1D("Valid_03", "Total Pz    ",100,-100.0, 100.0);
		hValid[ 4] = new TH1D("Valid_04", "Total P     ",100,   0.0, 100.0);
		hValid[ 5] = new TH1D("Valid_05", "phi of Track",100,  -3.2,   3.2);
		hValid[ 6] = new TH1D("Valid_06", "Vr  of Track",100, -10.0,  10.0); 
		hValid[ 7] = new TH1D("Valid_07", "Vz  of Track",100, -10.0,  10.0); 
		hValid[ 8] = new TH1D("Valid_08", "# of Tracks ",100,  -1.0,  99.0);
		hValid[ 9] = new TH1D("Valid_09", "# of PFOs   ",150,  -1.0, 149.0);
		hValid[10] = new TH1D("Valid_10", "Px of PFO   ",100, -50.0,  50.0);
		hValid[11] = new TH1D("Valid_11", "Py of PFO   ",100, -50.0,  50.0);
		hValid[12] = new TH1D("Valid_12", "Pz of PFO   ",100, -50.0,  50.0);
		hValid[13] = new TH1D("Valid_13", "E  of PFO   ",100,   0.0, 100.0);
		hValid[14] = new TH1D("Valid_14", "cosTheta PFO",100,  -1.0,   1.0);
		hValid[15] = new TH1D("Valid_15", "Inv. Mass ll",100,  20.0, 120.0);
		hValid[16] = new TH1D("Valid_16", "# of hits   ",100,   0.0, 200.0);
		hValid[17] = new TH1D("Valid_17", "Calo Energy ",100,   0.0, 100.0);
		hValid[18] = new TH1D("Valid_18", "cone Energy ",100,   0.0, 100.0);
		hValid[19] = new TH1D("Valid_19", "cone EN     ",100,   0.0, 100.0);
		hValid[20] = new TH1D("Valid_20", "InvMass All ",300,   0.0, 300.0);
		hValid[21] = new TH1D("Valid_21", "InvMass Hbb ",300,   0.0, 300.0);
		hValid[22] = new TH1D("Valid_22", "InvMass Hcc ",300,   0.0, 300.0);
		hValid[23] = new TH1D("Valid_23", "InvMass Hgg ",300,   0.0, 300.0);
		hValid[24] = new TH1D("Valid_24", "InvMass Hpp ",300,   0.0, 300.0);
		hValid[25] = new TH1D("Valid_25", "InvMass Huu ",300,   0.0, 300.0);
		hValid[26] = new TH1D("Valid_26", "InvMass Hzz ",300,   0.0, 300.0);
		hValid[27] = new TH1D("Valid_27", "InvMass Hww ",300,   0.0, 300.0);
		hValid[28] = new TH1D("Valid_28", "InvMass Htt ",300,   0.0, 300.0);
		//
		//
		hPando[ 0] = new TH1D("Pando_00", "Total Energy",100,   0.0, 300.0);
		hPando[ 1] = new TH1D("Pando_01", "Total Px    ",100,-100.0, 100.0);
		hPando[ 2] = new TH1D("Pando_02", "Total Py    ",100,-100.0, 100.0);
		hPando[ 3] = new TH1D("Pando_03", "Total Pz    ",100,-100.0, 100.0);
		hPando[ 4] = new TH1D("Pando_04", "Total P     ",100,   0.0, 100.0);
		hPando[ 5] = new TH1D("Pando_05", "phi of Track",100,  -3.2,   3.2);
		hPando[ 6] = new TH1D("Pando_06", "Vr  of Track",100, -10.0,  10.0); 
		hPando[ 7] = new TH1D("Pando_07", "Vz  of Track",100, -10.0,  10.0); 
		hPando[ 8] = new TH1D("Pando_08", "# of Tracks ",100,  -1.0,  99.0);
		hPando[ 9] = new TH1D("Pando_09", "# of PFOs   ",150,  -1.0, 149.0);
		hPando[10] = new TH1D("Pando_10", "Px of PFO   ",100, -50.0,  50.0);
		hPando[11] = new TH1D("Pando_11", "Py of PFO   ",100, -50.0,  50.0);
		hPando[12] = new TH1D("Pando_12", "Pz of PFO   ",100, -50.0,  50.0);
		hPando[13] = new TH1D("Pando_13", "E  of PFO   ",100,   0.0, 100.0);
		hPando[14] = new TH1D("Pando_14", "cosTheta PFO",100,  -1.0,   1.0);
		hPando[15] = new TH1D("Pando_15", "Inv. Mass ll",100,  20.0, 120.0);
		hPando[16] = new TH1D("Pando_16", "# of hits   ",100,   0.0, 200.0);
		hPando[17] = new TH1D("Pando_17", "Calo Energy ",100,   0.0, 100.0);
		hPando[18] = new TH1D("Pando_18", "cone Energy ",100,   0.0, 100.0);
		hPando[19] = new TH1D("Pando_19", "cone EN     ",100,   0.0, 100.0);
		hPando[20] = new TH1D("Pando_20", "InvMass All ",300,   0.0, 300.0);
		hPando[21] = new TH1D("Pando_21", "InvMass Hbb ",300,   0.0, 300.0);
		hPando[22] = new TH1D("Pando_22", "InvMass Hcc ",300,   0.0, 300.0);
		hPando[23] = new TH1D("Pando_23", "InvMass Hgg ",300,   0.0, 300.0);
		hPando[24] = new TH1D("Pando_24", "InvMass Hpp ",300,   0.0, 300.0);
		hPando[25] = new TH1D("Pando_25", "InvMass Huu ",300,   0.0, 300.0);
		hPando[26] = new TH1D("Pando_26", "InvMass Hzz ",300,   0.0, 300.0);
		hPando[27] = new TH1D("Pando_27", "InvMass Hww ",300,   0.0, 300.0);
		hPando[28] = new TH1D("Pando_28", "InvMass Htt ",300,   0.0, 300.0);

		hRes  [ 0] = new TH1D("hRes_ee",  "InvMass ee  ",1000,   0.0, 100.0);
		hRes  [ 1] = new TH1D("hRes_uu",  "InvMass uu  ",1000,   0.0, 100.0);
	}
	double _isrpzmax = 125;
	b = (double) 0.00464564*( std::log(m_ECM*m_ECM*3814714.)-1. );
	//= 2*alpha/pi*( ln(s/m_e^2)-1 )
	ISRPzMaxB = std::pow((double)_isrpzmax,b);
	//
}

void FSClasserProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void FSClasserProcessor::processEvent( LCEvent * evt ) { 
	_nEvt ++ ;
	_FD   =-1.;
	_Pmax = 0 ;
	_Emax = 0 ;
	if( m_debug>1 ) printf("\n");
	if( m_debug>1 || _nEvt%100==0 ) printf("Run & event =  %8d %8d\n",_nRun, _nEvt);
	m_HiggsBoostVector = TVector3(0,0,0); 
	//
	LCCollection    *colMC          = evt->getCollection(_colMCP );
	LCCollection    *col2Jets       = NULL; 
	LCCollection    *col4Jets       = NULL; 
	LCCollection    *colTaus        = NULL; 
	LCCollection    *colPFOs        = NULL; 
	LCCollection    *colIsoLeps     = NULL; 
	LCCollection    *colIsoZLeps    = NULL; 
	LCCollection    *colPFOPandora  = NULL; 
	LCCollection    *colMCTL        = NULL;
	if( m_full>0 && m_makeplots>0 )  colPFOPandora = evt->getCollection( "PandoraPFOs" );
	//	
	if( _colPFOs != "" ) colPFOs = evt->getCollection( _colPFOs );
	//	
	if( (m_findTaup ||  m_findTaum) && _colTaus != ""  ){
		colTaus  = evt->getCollection(_colTaus);
	}
	//	
	int nIsoLep = 0;
	if( _colIsoLeps != "") 
	{
		colIsoLeps  = evt->getCollection(_colIsoLeps);
		nIsoLep = colIsoLeps->getNumberOfElements();
		if( m_debug>1 ) { 
			printf("No. of elements in  LepCol = %3d\n", colIsoLeps->getNumberOfElements() );
		}
	}
	//	
	if( _colIsoZLeps != "") 
	{
		colIsoZLeps  = evt->getCollection(_colIsoZLeps);
		if( m_debug>1 ) { 
			printf("No. of elements in ZLepCol = %3d\n", colIsoZLeps->getNumberOfElements() );
		}
	}else if( _colIsoLeps != "") {

		colIsoZLeps  = colIsoLeps; 
	}
	//	
	if( _col2Jets != ""){
		 try{
			 col2Jets  = evt->getCollection(_col2Jets);
		 }catch(...){
			throw marlin::SkipEventException(this);
		 }
		 if( col2Jets ){
			if( m_debug>1 )  
				printf("No. of elements in  2JetCol = %3d\n", col2Jets->getNumberOfElements()    );
		}else{
			throw marlin::SkipEventException(this);
		}
	}
	//	
	if( _col4Jets != ""){
		 try{
			 col4Jets  = evt->getCollection(_col4Jets);
		 }catch(...){
			throw marlin::SkipEventException(this);
		 }
		 if( col4Jets ){
			if( m_debug>1 )  
				printf("No. of elements in  4JetCol = %3d\n", col4Jets->getNumberOfElements()    );
		}else{
			throw marlin::SkipEventException(this);
		}
	}
	//	
	if( m_matchmc>0 ) {
		colMCTL = evt->getCollection( _colMCTL );
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
	if(m_debug>1){ printf("nhfs = %3d\n", _nhfs);}
	//
	BuildFullSimParticleList( colMC, col2Jets, col4Jets, colMCTL,   colPFOs, colPFOPandora, colIsoZLeps, colTaus );
	/*
	   sort( photonList.begin(), photonList.end(), Sort_PFOs_Ed );
	   sort(     epList.begin(),     epList.end(), Sort_PFOs_Ed );
	   sort(     emList.begin(),     emList.end(), Sort_PFOs_Ed );
	   sort(    mupList.begin(),    mupList.end(), Sort_PFOs_Ed );
	   sort(    mumList.begin(),    mumList.end(), Sort_PFOs_Ed );
	   sort(    jetList.begin(),    jetList.end(), Sort_PFOs_Ed );
	*/
	if(m_debug>1){
		if (m_findJets   ) printf("No. of 2jets      = %3d\n",  (int)    jetList.size()) ;
		if (m_findJets   ) printf("No. of 4jets      = %3d\n",  (int)   jet4List.size()) ;
		if (m_findPhoton ) printf("No. of photons    = %3d\n",  (int) photonList.size()) ;
		if (m_findEp     ) printf("No. of electrons  = %3d\n",  (int)     epList.size()) ;
		if (m_findEm     ) printf("No. of positrons  = %3d\n",  (int)     emList.size()) ;
		if (m_findMup    ) printf("No. of muon+      = %3d\n",  (int)    mupList.size()) ;
		if (m_findMum    ) printf("No. of muon-      = %3d\n",  (int)    mumList.size()) ;
		if (m_findTaup   ) printf("No. of  tau+      = %3d\n",  (int)   taupList.size()) ;
		if (m_findTaum   ) printf("No. of  tau-      = %3d\n",  (int)   taumList.size()) ;
		if (1            ) printf("No. of PFOs       = %3d\n",        nPFOs            ) ;
	}
	//
	DictPList.insert(pair<string, vector<FSParticle*> > ("e+"     ,     epList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("e-"     ,     emList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("mu+"    ,    mupList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("mu-"    ,    mumList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("tau+"   ,   taupList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("tau-"   ,   taumList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("gamma"  , photonList )); 
	DictPList.insert(pair<string, vector<FSParticle*> > ("jet"    ,    jetList )); 
	// 
	//********************************************************************
	//
	//   LOOP OVER FINAL STATES
	//
	//********************************************************************
	//
	for (unsigned int ifs = 0; ifs < m_FSInfoVector.size(); ifs++){

		m_cutpass[ifs][0]++;

		// get information about this final state
		FSInfo* fsinfo                = m_FSInfoVector[ifs];
		//fsinfo->Print();
		NTupleHelper* NT              = fsinfo->NT();
		bool exclusive                = fsinfo->exclusive();
		bool inclusive                = fsinfo->inclusive();
		vector<string> particleNames  = fsinfo->particleNames();
		vector<int>    particleStatus = fsinfo->particleStatus();

		// check number of tracks
		if (exclusive){
			//if (evtRecEvent->totalCharged() > fsinfo->nChargedParticles()+m_maxExtraTracks) continue;
		}
		m_cutpass[ifs][1]++;

		//********************************************************************
		//
		//   PUT TOGETHER ALL COMBINATIONS OF PARTICLES FOR THIS FINAL STATE
		//
		//********************************************************************
		vector< vector< FSParticle* > > pCombos; pCombos.clear();
		for (unsigned int ifsp = 0; ifsp < particleNames.size(); ifsp++){

			int   status = particleStatus[ifsp];
			vector< FSParticle* >  pList ; pList.clear(); 
			if( status == 1 ) pList  = DictPList[particleNames[ifsp]];
			else{
				pList.push_back(new FSParticle(particleNames[ifsp], true));
				ParticleTrash.push_back(pList[0]); // save address for memory freeing 
			}
			vector< vector< FSParticle* > > pCombosTemp = pCombos; pCombos.clear();
			if  ( pList.size() == 0 ) break;

			for ( unsigned int ipl = 0; ipl < pList.size(); ipl++ ){

				if ( ifsp == 0 ) {

					vector< FSParticle* > combo1; combo1.clear();
					combo1.push_back(pList[ipl]); 
					if ( checkCombination( combo1, (combo1.size()==particleNames.size()), fsinfo->inclusive()||fsinfo->missingMassFit() ) )
						pCombos.push_back(combo1);

				} else {

					for (unsigned int itc = 0; itc < pCombosTemp.size(); itc++){
						vector< FSParticle* > combo2 = pCombosTemp[itc];
						bool duplicate = false;
						if  ( !pList[ipl]->missed() ) {  
							for (unsigned int ic = 0; ic < combo2.size(); ic++){
								if ( combo2[ic]->missed() ) continue;
								bool shared = pList[ipl]->duplicate(combo2[ic], m_full);
								if ( shared ) { duplicate = true; break; }
							}
							if(m_debug>4){
								printf("duplicated one ? %2d\n", duplicate );
							}
							if ( duplicate ) continue;
						}
						combo2.push_back(pList[ipl]); 
						if ( checkCombination( combo2, (combo2.size()==particleNames.size()), fsinfo->inclusive()||fsinfo->missingMassFit() ) )
							pCombos.push_back(combo2);

						if ( pCombos.size()>100000 ) break;
					}

				}
				if (pCombos.size() > 100000) break;

			}

			if (pCombos.size() > 100000) break;

		}
		m_cutpass[ifs][2]++;
		//	
		if(m_debug > 1)
			printf("1: Channel: %20s, No. = %10d, Size of pCombos is %4d\n",
					fsinfo->FSName().data(), ifs, (int)pCombos.size()); 
		if( pCombos.size() ==  0) continue;
		if( pCombos.size() >99999) {
			printf("FSClasser WAENING: possible combination is too many, %5d, Skipping ... \n",(int)pCombos.size());
			cout << "\tRun      = " << _nRun;
			cout << ", Event    = " << _nEvt<< endl;
			cout << "\t\tFINAL STATE = " << fsinfo->FSName() << endl;
			continue;
		}
		m_cutpass[ifs][3]++;
		//********************************************************************
		//
		//   LOOP OVER PARTICLE COMBINATIONS
		//
		//********************************************************************
		//********************************************************************
		//
		for (unsigned int i = 0; i < pCombos.size(); i++){
			m_cutpass[ifs][4]++;
			if(m_debug > 2)
				printf("cutpass[%2d][4] =  %5d %7d\n",i,(int)m_cutpass[ifs][4], _nEvt);
			vector<FSParticle*> combo = pCombos[i];
			//********************************************************************
			//
			//   CUT ON THE TOTAL ENERGY AND MOMENTUM TO SAVE TIME
			//
			//********************************************************************

			totalP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
			for (unsigned int t = 0; t < combo.size(); t++){
				if(!combo[t]->missed() ) totalP4 += combo[t]->rawFourMomentum();
			}
			for (unsigned int t = 0; t < combo.size(); t++){
				if( combo[t]->missed() ) combo[t]->setRawFourMomentum (m_ecms-totalP4);
			}
			double missingEnergy = (totalP4-m_ecms). E();
			double missingMass   = (totalP4-m_ecms). M();
			double missingMass2  = (totalP4-m_ecms).M2();
			if(m_debug > 2) {
				printf("3: VisP4 = %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",pVis   .Px(), pVis   .Py(), pVis   .Pz(), pVis   .E(), m_mcTruthHelper->MCTotalEnergy(), m_mcTruthHelper->MCMissingEnergy() );
				printf("3: TotP4 = %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",totalP4.Px(), totalP4.Py(), totalP4.Pz(), totalP4.E(), m_mcTruthHelper->MCTotalEnergy(), m_mcTruthHelper->MCMissingEnergy() );
				printf("4: miss  = %10.4f %10.4f %10.4f\n",fabs(missingEnergy), missingMass, fsinfo->missingMassValue());
			}
			//
			if (exclusive ){
				if( fsinfo->missingMassFit() ){
					if(fabs(missingMass-fsinfo->missingMassValue())> 1000)continue;
				}	
				else if (fabs(missingEnergy)                      >10000)continue;
			}
			m_cutpass[ifs][5]++;
			//
			if (exclusive && !fsinfo->missingMassFit() ){
				double missingP = (m_ecms-totalP4).Rho();	
				if (fabs(missingP)>10000) continue;
			}
			m_cutpass[ifs][6]++;
			//
			if (exclusive){
				if ( fabs(missingMass2) > 10000) continue;
			}
			if (inclusive){
				//if (missingMass2 < m_lowerMissingMass2) continue;
				//if (missingMass2 > m_upperMissingMass2) continue;
			}
			m_cutpass[ifs][7]++;

			//********************************************************************
			//
			//   MAKE CUTS ON INTERMEDIATE PARTICLE COMBINATIONS (RAW 4-VECTORS)
			//
			//********************************************************************

			if (!(fsinfo->evaluateFSCuts(combo, m_ecms,"Raw"))) continue;
			m_cutpass[ifs][8]++;
			//********************************************************************
			//
			//   DO THE KINEMATIC FITTING
			//
			//********************************************************************
			double prob    = 0.0;
			double chi2    = 99999;
			double chi2ndf = 99999;
			double nitr    = -1;
			vector<double> scale; scale.clear();

			if ( exclusive && m_kmfit <0 && combo.size()==4){
				////////////////////////////////////////////////////////////////
				// 2 fermions and 2 photons energy and momentum balance using matrix method
				// A * X = B,  X=(a,b,c,d) to be determined, B=(sqrtS,0,0,0),
				// A is matrix of 4-momentum of 4 jets
				//          A            *   X   =      B
				//  { E1  E2  E3  E4  }    { a }    { sqrtS }
				//  | Px1 Px2 Px3 Px4 |  * | b | =  |   0   |
				//  | Py1 Py2 Py3 Py4 |    | c |    |   0   |
				//  { Pz1 Pz2 Pz3 Pz4 }    { d }    {   0   }
				//
				// then X = invert(A) * B
				// The idea is to find a set of parameters (a,b,c,d) to reweight the Energy, Px, Py, Pz of 4 jets
				// a*E1  + b*E2  + c*E3  + d*E4  = sqrtS (total e = 250)
				// a*Px1 + b*Px2 + c*Px3 + d*Px4 =  0    (total px = 0)
				// a*Py1 + b*Py2 + c*Py3 + d*Py4 =  0    (total py = 0)
				// a*Pz1 + b*Pz2 + c*Pz3 + d*Pz4 =  0    (total pz = 0)
				// ----------------------------------------------------------
				// discussion with Manqi Ruan (IHEP) and Jianping Dai (SJTU)
				// written by Haijun Yang @ Shanghai Jiao Tong University, 2014/06/17
				// ----------------------------------------------------------
				unsigned int msize = combo.size();
				TMatrixD H0 = THilbertMatrixD(msize,msize);

				for (unsigned int j = 0; j < combo.size(); j++){
					FSParticle* fsp = combo[j];
					H0[0][j] = fsp->rawFourMomentum().E();
					H0[1][j] = fsp->rawFourMomentum().X();
					H0[2][j] = fsp->rawFourMomentum().Y();
					H0[3][j] = fsp->rawFourMomentum().Z();
				}
				double det1;
				TMatrixD H1 = H0;
				TMatrixD H2 = H1.Invert(&det1);
				// H1.Print();
				// H2.Print();
				// Get the maximum off-diagonal matrix value . One way to do this is to set diagonal to zero .
				// TMatrixD U1(H1,TMatrixD::kMult,H_square);
				// TMatrixDDiag diag1(U1); diag1 = 0.0;
				// const Double_t U1_max_offdiag = (U1.Abs()).Max();
				// cout << "  Maximum off-diagonal = " << U1_max_offdiag << endl;
				// cout << "  Determinant          = " << det1 <<endl;
				for (unsigned int j = 0; j < combo.size(); j++){
					FSParticle* fsp = combo[j];
					double fs = H2[j][0]*m_ECM; 
					scale.push_back( fs );
					fsp->setFitFourMomentum( fsp->rawFourMomentum()*fs );
				}
				//double cut = 0.5;
				//if ( fabs(a-1.)>cut || fabs(b-1.)>cut || fabs(c-1.)>cut || fabs(d-1.)>cut ) { a = 1.; b = 1.; c = 1.; d = 1.; }
				//cout << " a/b/c/d = " << a << " " << b << " " << c << " " << d << endl;
			}
			// initialize
			if ( exclusive && m_kmfit >0 ) {
				//
				vector<JetFitObject>  fitjets;
				for (unsigned int j = 0; j < combo.size(); j++) {
					if ( combo[j]->missed() ){}
					else{
						fitjets.push_back( *(combo[j]->jetfitobject()) );
					} 
				}
				//
				PConstraint pxc (1, 0, 0, 0);
				PConstraint pyc (0, 1, 0, 0);
				PConstraint pzc (0, 0, 1, 0);
				PConstraint pec (0, 0, 0, 1, m_ECM);
				if(0){
					ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
					pxc.addToFOList (*(photon));
					pyc.addToFOList (*(photon));
					pzc.addToFOList (*(photon));
					pec.addToFOList (*(photon));
				}
				//
				NewFitterGSL kmfit;	
				// OPALFitterGSL kmfit;	
				// add final state particles
				int kinFitTrack = 0;
				int kinFitDOF   = 0;
				int kinFitConstraint = 0;
				//vector<int> exList; exList.clear();
				for (unsigned int j = 0; j < combo.size(); j++){
					FSParticle* fsp = combo[j];
					if ( fsp->missed() ){
					} else{
						/*
						   kmfit.addFitObject( *( combo[j]->jetfitobject()) ); 
						   pxc.addToFOList   ( *( combo[j]->jetfitobject()) );
						   pyc.addToFOList   ( *( combo[j]->jetfitobject()) );
						   pzc.addToFOList   ( *( combo[j]->jetfitobject()) );
						   pec.addToFOList   ( *( combo[j]->jetfitobject()) );
						   */
						//	
						kmfit.addFitObject(  fitjets[j]  ); 
						pxc.addToFOList   (  fitjets[j]  );
						pyc.addToFOList   (  fitjets[j]  );
						pzc.addToFOList   (  fitjets[j]  );
						pec.addToFOList   (  fitjets[j]  );
						//
						kinFitTrack++;
					}
				}
				if ( fsinfo->constrain4Mom() ){
					kmfit.addConstraint (pxc);
					kmfit.addConstraint (pyc);
					kmfit.addConstraint (pzc);
					kmfit.addConstraint (pec);
					kinFitDOF += 4;
					kinFitConstraint++;
				}
				// do the fit and make a very loose cut on the resulting chi2
				prob    = kmfit.fit();
				chi2    = kmfit.getChi2();
				chi2ndf = kmfit.getChi2()/kinFitDOF;
				nitr    = kmfit.getIterations();
				//
				if (kinFitConstraint > 0){ 
					kinFitTrack = 0;
					for (unsigned int j = 0; j < combo.size(); j++){
						combo[j]->setFitFourMomentum(
								TLorentzVector(
									fitjets[j].getPx(), 
									fitjets[j].getPy(), 
									fitjets[j].getPz(), 
									fitjets[j].getE() )
								);
					}
					//if( m_debug>1) cout<<"E of ISR photon = "<<photon->getE()<<endl;;

				}
				fitjets.clear();
				//delete photon; 
				//********************************************************************
				//
				//   MAKE CUTS ON INTERMEDIATE PARTICLE COMBINATIONS (FIT 4-VECTORS)
				//
				//********************************************************************

				if (exclusive && !(fsinfo->evaluateFSCuts(combo, m_ecms,"Fit"))) continue;
			}
			m_cutpass[ifs][9]++;
			//********************************************************************
			//
			//   RECORD INFORMATION
			//
			//********************************************************************
			// record event level information

			if( m_debug>3 ) printf("2: F: fill  ... %3d %3d\n",i, ifs); 
			NT->fillDouble("Run",         _nRun);
			NT->fillDouble("Event",       _nEvt);
			NT->fillDouble("Weight",     m_wate);
			NT->fillDouble("ntrks",       ntrks);
			NT->fillDouble("nclus",       nclus);
			NT->fillDouble("nPFOs",       nPFOs);
			NT->fillDouble("ntrks_Pandora",ntrks_Pandora);
			NT->fillDouble("nclus_Pandora",nclus_Pandora);
			NT->fillDouble("nPFOs_Pandora",nPFOs_Pandora);
			NT->fillDouble("Pmax",        _Pmax);
			NT->fillDouble("Emax",        _Emax);
			NT->fillDouble("njets",       numberJets);
			NT->fillDouble("ntaus",       numberTaus);
			NT->fillDouble("nElec",       epList.size() + emList.size());
			NT->fillDouble("nMuon",       mupList.size()+mumList.size());
			NT->fillDouble("nIsoLep",     nIsoLep);
			NT->fillDouble("nGamma",      photonList.size());
			NT->fillDouble("VisEn",       pVis.E());
			NT->fillDouble("VisPx",       pVis.X());
			NT->fillDouble("VisPy",       pVis.Y());
			NT->fillDouble("VisPz",       pVis.Z());
			NT->fillDouble("VisMass",     pVis.M());
			//	
			//
			NT->fillDouble("y12",   _y12  );
			NT->fillDouble("y23",   _y23  );
			NT->fillDouble("y34",   _y34  );
			NT->fillDouble("y45",   _y45  );
			NT->fillDouble("y56",   _y56  );
			NT->fillDouble("y67",   _y67  );
			if( m_EventShape>0){
				if( m_boostFS.size()>=4 && (m_boostFS.find("_") != string::npos) ) 
					m_HiggsBoostVector = fsinfo->HiggsBoostVector(combo, m_boostFS);
				//m_HiggsBoostVector.Print();
				if( m_LinearSphericity>0) Lsphericity();
				else                       sphericity();
				thrust();
				NT->fillDouble("Sphericity", m_Sphericity  );
				NT->fillDouble("Aplanarity", m_Aplanarity  );
				NT->fillDouble("C",          m_C           );
				NT->fillDouble("D",          m_D           );
				NT->fillDouble("Thrust",     m_Thrust      );
				NT->fillDouble("Major",      m_Major       );
				NT->fillDouble("Minor",      m_Minor       );
				NT->fillDouble("ThrustTheta",cos(m_Theta)  );
				NT->fillDouble("ThrustPhi",  m_Phi         );
				NT->fillDouble("ThrustEDM",  m_ThrEDM      );
				NT->fillDouble("MajorTheta", cos(m_Thetap) );
				NT->fillDouble("MajorPhi",   m_Phip        );
				NT->fillDouble("MajorEDM",   m_ThrEDMp     );
				NT->fillDouble("FD",         FDevent()     );
			}
			
			
			//if(m_savemc)NT->fillMCTruth(m_mcTruthHelper);
			if(m_luxury && m_full>0){
				NT->fillTracks(raw_trackvec);
			}
			if(m_luxury){
				NT->fillMCTruth(m_mcTruthHelper);
				NT->fillPFOs(raw_PFOs, colMCTL);
			}
			//
			NT->fillDouble("nhfs" ,    _nhfs  );
			NT->fillDouble("VisEnMC",  m_mcTruthHelper->MCTotalEnergy()   );
			NT->fillDouble("MisEnMC",  m_mcTruthHelper->MCMissingEnergy() );

			if( m_debug>3 ) printf("2: F: fill  ...   global\n");

			// record odds and ends that have no other place
			NT->fillDouble("MissingMass2",   missingMass2 );
			NT->fillDouble("TotalP",         totalP4.Rho());
			//
			if (m_kmfit!=0 && (exclusive || fsinfo->constrain4Mom()) ){
				NT->fillDouble("Chi2",   chi2);
				NT->fillDouble("prob",   prob);
				NT->fillDouble("nitr",   nitr);
				NT->fillDouble("Chi2DOF",chi2ndf);
			}

			if (inclusive){
				NT->fillDouble("TotalEnergy",    totalP4.E()  );
				NT->fillDouble("TotalPx",        totalP4.Px() );
				NT->fillDouble("TotalPy",        totalP4.Py() );
				NT->fillDouble("TotalPz",        totalP4.Pz() );
			}
			if ( m_kmfit<0 && scale.size()>0) {
				NT->fillArray ("kScale"    , "idx_kmf", scale, scale.size());
			}
			//record particle level information
			vector<TLorentzVector> rawp4list, kmfp4list;
			rawp4list.clear();  
			kmfp4list.clear(); 
			vector<FSParticle*>    gammaPFOs, muonPFOs, elecPFOs;
			gammaPFOs.clear();
			muonPFOs .clear(); 
			elecPFOs .clear(); 
			//
			TLorentzVector rawp4gamma(0,0,0,0);
			for (unsigned int j = 0; j < photonList.size(); j++){
				rawp4gamma += photonList[j]->rawFourMomentum();
				int m=0;
				for (unsigned int mn = 0; mn < combo.size(); mn++){
					FSParticle* fsp = combo[mn];
					if ( fsp == photonList[j] ) m++;
				}
				if( m==0 ){
					gammaPFOs.push_back(photonList[j]);
				}
			}
			//
			for (unsigned int j = 0; j < mupList.size(); j++){
				int m=0;
				for (unsigned int mn = 0; mn < combo.size(); mn++){
					FSParticle* fsp = combo[mn];
					if ( fsp == mupList[j] ) m++;
				}
				if( m==0 ) {
					muonPFOs.push_back(mupList[j]);
				}
			}
			for (unsigned int j = 0; j < mumList.size(); j++){
				int m=0;
				for (unsigned int mn = 0; mn < combo.size(); mn++){
					FSParticle* fsp = combo[mn];
					if ( fsp == mumList[j] ) m++;
				}
				if( m==0 ) {
					muonPFOs.push_back(mumList[j]);
				}
			}
			//
			for (unsigned int j = 0; j < epList.size(); j++){
				int m=0;
				for (unsigned int mn = 0; mn < combo.size(); mn++){
					FSParticle* fsp = combo[mn];
					if ( fsp == epList[j] ) m++;
				}
				if( m==0 ) {
					elecPFOs.push_back(epList[j]);
				}
			}
			for (unsigned int j = 0; j < emList.size(); j++){
				int m=0;
				for (unsigned int mn = 0; mn < combo.size(); mn++){
					FSParticle* fsp = combo[mn];
					if ( fsp == emList[j] ) m++;
				}
				if( m==0 ) {
					elecPFOs.push_back(emList[j]);
				}
			}
			//
			char index[30];
			for (unsigned int j = 0; j < combo.size(); j++){
				FSParticle* fsp = combo[j];
				rawp4list.push_back(fsp->rawFourMomentum());
				if( m_kmfit!=0 && exclusive ) kmfp4list.push_back(fsp->fitFourMomentum());
				if( fsp->type()!=4)
					NT->fillPFO(fsp->pfo(), fsp, j+1,"Pfo", m_full );
				else            
					NT->fillJet(fsp            , j+1,"Jet", m_full , m_TagFlavor, m_kappa);
			}
			for (unsigned int j = 0; j < jet4List.size(); j++){
				FSParticle* fsp = jet4List[j];  
				NT->fillJet(fsp            , j+1,"AltJet", m_full , m_TagFlavor, m_kappa);
			}
	
			//
			NT->fill4Momentum("idx_raw","raw_", rawp4list, rawp4list.size());
			if(m_luxury){
				NT->fillDouble("RMassAllPhotons",   rawp4gamma.M() );
				NT->fillExtras(  elecPFOs,  "ExElec_",  11);
				NT->fillExtras(  muonPFOs,  "ExMuon_",  13);
				NT->fillExtras( gammaPFOs, "ExGamma_",  22);
			}
			if(exclusive && m_kmfit!=0)NT->fill4Momentum("idx_kmf","kmf_", kmfp4list, kmfp4list.size());
			if( rawp4list.size()>=1 && inclusive){
				for(unsigned int ki=0; ki<rawp4list.size(); ki++){
					sprintf(index ,"Rreco%d",ki+1);
					NT->fillDouble((string)index,(m_ecms-rawp4list[ki]).M2());
					if (exclusive &&  m_kmfit!=0 && kmfp4list.size()==rawp4list.size()) {
						sprintf(index ,"Kreco%d",ki+1);
						NT->fillDouble((string)index,(m_ecms-kmfp4list[ki]).M2());
					}
				}
			}
			if( rawp4list.size()>=2 ){
				for(unsigned int ki=0; ki<rawp4list.size()-1; ki++){
					for(unsigned int kj=ki+1; kj<rawp4list.size(); kj++){
						sprintf(index ,"RMass%d%d",ki+1,kj+1);
						NT->fillDouble((string)index,(rawp4list[ki]+rawp4list[kj]).M());
						sprintf(index ,"RCosTheta%d%d",ki+1,kj+1);
						NT->fillDouble((string)index,cos(rawp4list[ki].Angle(rawp4list[kj].Vect())));
						if (exclusive &&   m_kmfit!=0  ) {
							sprintf(index ,"KMass%d%d",ki+1,kj+1);
							NT->fillDouble((string)index,(kmfp4list[ki]+kmfp4list[kj]).M());
							sprintf(index ,"KCosTheta%d%d",ki+1,kj+1);
							NT->fillDouble((string)index,cos(kmfp4list[ki].Angle(kmfp4list[kj].Vect())));
						}
						sprintf(index ,"Rreco%d%d",ki+1,kj+1);
						NT->fillDouble((string)index,(m_ecms-rawp4list[ki]-rawp4list[kj]).M());
						if (exclusive &&  m_kmfit!=0 && kmfp4list.size()==rawp4list.size()) {
							sprintf(index ,"Kreco%d%d",ki+1,kj+1);
							NT->fillDouble((string)index,(m_ecms-kmfp4list[ki]-kmfp4list[kj]).M2());
						}
					}
				}
			}
			if( rawp4list.size()>=3 ){
				for(unsigned int ki=0; ki<rawp4list.size()-2; ki++){
					for(unsigned int kj=ki+1; kj<rawp4list.size()-1; kj++){
						for(unsigned int kk=kj+1; kk<rawp4list.size(); kk++){
							sprintf(index ,"RMass%d%d%d",ki+1,kj+1,kk+1);
							NT->fillDouble((string)index,(rawp4list[ki]+rawp4list[kj]+rawp4list[kk]).M());
							if (exclusive &&  m_kmfit!=0  ) {
								sprintf(index ,"KMass%d%d%d",ki+1,kj+1,kk+1);
								NT->fillDouble((string)index,(kmfp4list[ki]+kmfp4list[kj]+kmfp4list[kk]).M());
							}
							sprintf(index ,"Rreco%d%d%d",ki+1,kj+1,kk+1);
							NT->fillDouble((string)index,(m_ecms-rawp4list[ki]-rawp4list[kj]-rawp4list[kk]).M());
							if (exclusive &&  m_kmfit!=0 && kmfp4list.size()==rawp4list.size()) {
								sprintf(index ,"Kreco%d%d%d",ki+1,kj+1,kk+1);
								NT->fillDouble((string)index,(m_ecms-kmfp4list[ki]-kmfp4list[kj]-kmfp4list[kk]).M2());
							}
						}
					}
				}
			}
			if( rawp4list.size()>=4 ){
				for(unsigned int kl=0; kl<rawp4list.size()-3; kl++){
					for(unsigned int ki=kl+1; ki<rawp4list.size()-2; ki++){
						for(unsigned int kj=ki+1; kj<rawp4list.size()-1; kj++){
							for(unsigned int kk=kj+1; kk<rawp4list.size(); kk++){
								sprintf(index ,"RMass%d%d%d%d",kl+1,ki+1,kj+1,kk+1);
								NT->fillDouble((string)index,(rawp4list[kl]+rawp4list[ki]+rawp4list[kj]+rawp4list[kk]).M());
								if (exclusive &&   m_kmfit!=0  ) {
									sprintf(index ,"KMass%d%d%d%d",kl+1,ki+1,kj+1,kk+1);
									NT->fillDouble((string)index,(kmfp4list[kl]+kmfp4list[ki]+kmfp4list[kj]+kmfp4list[kk]).M());
								}
								sprintf(index ,"Rreco%d%d%d%d",kl+1,ki+1,kj+1,kk+1);
								NT->fillDouble((string)index,(m_ecms-rawp4list[kl]-rawp4list[ki]-rawp4list[kj]-rawp4list[kk]).M());
								if (exclusive &&  m_kmfit!=0 && kmfp4list.size()==rawp4list.size()) {
									sprintf(index ,"Kreco%d%d%d%d",kl+1,ki+1,kj+1,kk+1);
									NT->fillDouble((string)index,(m_ecms-kmfp4list[kl]-kmfp4list[ki]-kmfp4list[kj]-kmfp4list[kk]).M2());
								}
							}
						}
					}
				}
			}

			if ( fsinfo->hasParticle("tau+") &&  fsinfo->hasParticle("tau-") ){
				vector< vector<unsigned int> > indices = fsinfo->submodeIndices("0_0110000");
				double missingx = -totalP4.X();
				double missingy = -totalP4.Y();
				for (unsigned int i = 0; i < indices.size(); i++){
					vector<unsigned int> indexCombo = indices[i];
					if( indexCombo.size()==2){
						TLorentzVector p4plus  = combo[indexCombo[0]]->rawFourMomentum();
						TLorentzVector p4minus = combo[indexCombo[1]]->rawFourMomentum();
						double MassTT = getMassColApp( p4plus, p4minus, missingx, missingy  );
						sprintf(index ,"MassTT_%1d",i+1);
						NT->fillDouble((string)index,MassTT);
					}
				}
			}	
			// write the tree
			m_CutPass[ifs]++;
			NT->write();
			//
			FreeAll(rawp4list);
			FreeAll(kmfp4list);
			FreeAll(gammaPFOs);
			FreeAll( muonPFOs); 
			FreeAll( elecPFOs);
		}
	}
	//
	CleanEvt();
}

void FSClasserProcessor::check( LCEvent * evt ) { 

}

void FSClasserProcessor::MakePlots(){

	if ( m_savehis ){
		for(int i=0; i<=28; i++) {
			if(hValid[i] && hValid[i]->Integral()>0){
				TH1D *h1 = new TH1D("h1","",hValid[i]->GetNbinsX(), hValid[i]->GetXaxis()->GetXmin(), hValid[i]->GetXaxis()->GetXmax());
				TH1D *h2 = 0; 
				char filename[256], title[256];
				sprintf(filename,"figs/%s", (hValid[i]->GetName()));
				sprintf(title,"%s", (hValid[i]->GetTitle()));
				//
				PlotDataMC(filename, 
						h1       , (char*)"",
						hValid[i], (char*)"Arbor",
						hPando[i], (char*)"Pandora",
						h2       , (char*)"",
						h2       , (char*)"",
						true, false, false, title
					  );
				delete h1;
			}
		}
	}

}

void FSClasserProcessor::end(){
	if ( m_savehis ){
		if(m_full>0 && m_makeplots>0)MakePlots();
		for(int i=0; i<=28; i++) {
			if(hValid[i]){
				hValid[i]->Write();
				delete hValid[i];

				hPando[i]->Write();
				delete hPando[i];
			}

		}
		delete hRes[0];
		delete hRes[1];
		delete m_ntVal;
	}
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


bool FSClasserProcessor::IsPhoton( ReconstructedParticle *par, LCCollection *pfoCol ) {

	const double cosCone  = TMath::Cos(0.2); //11.5 grads

	double ecalE      = 0.;
	double hcalE      = 0.;
	int nConeCharged  = 0 ;
	int nConeNeutral  = 0 ;
	/* get deposited energies in eCal and hCal*/
	double cale[2]={0.0, 0.0};
	getCalEnergy(par, cale);
	double coneEC     = cale[0];
	double coneEN     = cale[1];
	double totCalE   = ecalE + hcalE;
	/* get number (and their energies) of charged/neutral particles inside the cone */
	getConeNeutraChargeInfo(par, pfoCol, cosCone, nConeCharged, nConeNeutral, coneEC, coneEN);
	//double coneE     = coneEC + coneEN;
	//                                                                                                                                                                                                   
	TVector3 momentum   = par->getMomentum();
	double pmag         = momentum.Mag();
	if(m_debug>2) printf("photon energies %8.2f %8.2f %8.2f\n",ecalE,hcalE, pmag);
	//double pt           = getTransverseMomentum(par);
	Double_t fEpsilon   = 1.E-10;
	/* true if the par is a photon from Higgs -> gamma gamma */
	bool neutralPar            = (par->getCharge() == 0);
	// asking isolation before correct for radiation decrease efficiency
	//bool isolated              = (coneE   < 5); // 5 is quite loose
	//bool highPt                = (pt      > 25);
	bool highE                 = (totCalE > 5);
	bool ecalE_over_Etot       = (ecalE/(ecalE+hcalE+fEpsilon)  > 0.9 );
	bool totCalE_over_pmag_min = ((ecalE+hcalE)/pmag            > 0.8 );
	bool totCalE_over_pmag_max = ((ecalE+hcalE)/pmag            < 1.2 );
	bool caloCuts              = ( ecalE_over_Etot && totCalE_over_pmag_min && totCalE_over_pmag_max);
	//        
	bool isPhoton = false;
	if(m_debug>2) printf("photon cuts %2d %2d %2d\n",neutralPar,highE,caloCuts);
	isPhoton = (  
			neutralPar     &&
			highE          &&  // energetic photons      
			caloCuts       );
	//
	return isPhoton;
}

void FSClasserProcessor::getCalEnergy(ReconstructedParticle *pfo, double *cale)
{
	double ecal = 0;
	double hcal = 0;
	std::vector<lcio::Cluster*> clusters = pfo->getClusters();
	
	if(m_debug>3) printf("%3d clusters in this pfo \n", (int) clusters.size() );

	for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
			iCluster!=clusters.end();
			++iCluster) {
		ecal += (*iCluster)->getSubdetectorEnergies()[0];
		hcal += (*iCluster)->getSubdetectorEnergies()[1];
		ecal += (*iCluster)->getSubdetectorEnergies()[3];
		hcal += (*iCluster)->getSubdetectorEnergies()[4];
		printf("calo  energies %7.3f %7.3f \n", ecal, hcal );
	}
	cale[0] = ecal;
	cale[1] = hcal;

	if(m_debug>3) printf("h/e calo  energies %7.3f %7.3f \n", ecal, hcal );

}


void FSClasserProcessor::getConeNeutraChargeInfo(ReconstructedParticle *part, LCCollection *collPfo, 
		double cosCone, int &nConeCharged, int &nConeNeutral,
		double &coneEC, double &coneEN)
/*
   Parameter are Reconstructed particle who defne the cone
   PandoraPfo collection
   consine of the angle defining the size of the cone
   and the variables we want to get: 
   ===
   number of charged/neutral particles inside cone
   and their energies
   */
{

	/*get cone information*/
	std::vector<lcio::ReconstructedParticle*> conePFOs;
	//Double_t coneEnergy = getConeEnergy(part,collPfo,cosCone,conePFOs);
	//int nConePFOs     = conePFOs.size(); /* number of particles inside cone */

	// check that they are initialized to zero
	if (nConeCharged != 0 || nConeNeutral != 0 || coneEC != 0. || coneEN != 0.) {
		std::cerr << " >>> WARNING: following variables should be zero, but they are not:" << std::endl
			<< " nConeCharged =  " << nConeCharged
			<< " nConeNeutral = "  << nConeNeutral
			<< " coneEC = "        << coneEC
			<< " coneEN = "        << coneEN << std::endl;

		std::cout << " Seting them to zero ... " << std::endl;
		nConeCharged = 0;
		nConeNeutral = 0;
		coneEC       = 0.;
		coneEN       = 0.;
	}



	for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=conePFOs.begin();iObj<conePFOs.end();++iObj) {
		if ((*iObj)->getCharge() != 0) {
			nConeCharged += 1;
			coneEC += (*iObj)->getEnergy();

		} else {
			nConeNeutral += 1;
			coneEN += (*iObj)->getEnergy();
			//std::cout << " energy: " << (*iObj)->getEnergy() 
			//          << " px: " << getMomentumX((*iObj)) << std::endl;

		}
	}

}

double FSClasserProcessor::getTransverseMomentum(ReconstructedParticle *part)
{

	TVector3 momentum        = TVector3(part->getMomentum());
	//double momentumMagnitude = momentum.Mag();
	//double sinTheta          = TMath::Sqrt( 1 - momentum.CosTheta() * momentum.CosTheta() );
	//double pT                = momentumMagnitude * sinTheta;

	return momentum.Pt();
}

void
FSClasserProcessor::CleanEvt(){
	//********************************************************************
	//
	//         CLEAN UP MEMORY
	//
	//********************************************************************
	FreeDelAll( epList        );
	FreeDelAll( emList        );
	FreeDelAll( mupList       );
	FreeDelAll( mumList       );
	FreeDelAll( taupList      );
	FreeDelAll( taumList      );
	FreeDelAll( photonList    );
	FreeDelAll(  jetList      );
	FreeDelAll( jet4List      );
	FreeDelAll( ParticleTrash );
	//
	FreeAll   ( DictPList     );
	//
	delete m_mcTruthHelper; m_mcTruthHelper=NULL;
}



bool
FSClasserProcessor::checkCombination(const vector<FSParticle*>& combo, bool complete, bool inclusive){

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
//
//
void    
FSClasserProcessor::BuildFullSimParticleList(
		LCCollection    *colMC  , LCCollection    *col2Jets, LCCollection    *col4Jets, 
		LCCollection   *colMCTL , LCCollection    *colPFO , 
		LCCollection *colPFOPandora, LCCollection  *colIsoLeps,
		LCCollection  *colTaus 
		){
	Double_t fCosConeCut   = 0.98;        // the angle of cone around the direction of pfo
	//
	//Double_t fElectronCut1 = 0.5;         // lower edge of totalCalEnergy/momentum
	//Double_t fElectronCut2 = 1.3;         // upper edge of totalCalEnergy/momentum
	//Double_t fElectronCut3 = 0.9;         // lower edge of ecalEnergy/totalCalEnergy
	//Double_t fMuonCut1     = 0.3;         // upper edge of totalCalEnergy/momentum
	////Double_t fMuonCut3     = 1.2;       // lower edge of yoke energy
	//Double_t fMuonCut3     = -1.2;        // lower edge of yoke energy, yoke energy is 0,  why? needs check
	//Double_t fEpsilon      = 1.E-10;
	//Double_t c0Electron    = 12.2;        // used in llHH mode
	//Double_t c1Electron    = 0.87; 
	//Double_t c0Muon        = 12.6; 
	//Double_t c1Muon        = 4.62;
	vector<MCParticleImpl*>  partonList = m_mcTruthHelper->partonList();
	//
	m_Sphericity = -1, m_Aplanarity = -1, m_C = -1, m_D = -1, m_Thrust=-1, m_Theta=999, m_Phi=999, m_ThrEDM=9e9;  
	ntrks         = 0;
	nclus         = 0;
	nPFOs         = 0;
	ntrks_Pandora = 0;
	nclus_Pandora = 0;
	nPFOs_Pandora = 0;
	numberJets    = 0; 
	pVis          = TLorentzVector(0.,0.,0.,0.);
	pVis_Pandora  = TLorentzVector(0.,0.,0.,0.);
	//
	vector<double> pd_phi     , ab_phi     ;
	vector<double> pd_Vr      , ab_Vr      ;
	vector<double> pd_Vz      , ab_Vz      ;
	vector<double> pd_costheta, ab_costheta;
	vector<double> pd_energy  , ab_energy  ;
	vector<double> pd_nHits   , ab_nHits   ;
	vector<double> pd_calEn   , ab_calEn   ;
	vector<double> pd_ecalE   , ab_ecalE   ;
	vector<double> pd_hcalE   , ab_hcalE   ;
	vector<TLorentzVector>  track_pd_p4list, track_ab_p4list, neutral_pd_p4list, neutral_ab_p4list;
	//
	TLorentzVector ab_pvis_trk = TLorentzVector(0.,0.,0.,0.);
	TLorentzVector ab_pvis_neu = TLorentzVector(0.,0.,0.,0.);
	TLorentzVector pd_pvis_trk = TLorentzVector(0.,0.,0.,0.);
	TLorentzVector pd_pvis_neu = TLorentzVector(0.,0.,0.,0.);
	//
	if( colPFOPandora ){
		nPFOs_Pandora = colPFOPandora->getNumberOfElements();
		for (int i=0;i<nPFOs_Pandora;i++) {
			ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFOPandora->getElementAt(i));
			TrackVec   tckvec = recPart->getTracks();
			ClusterVec cluvec = recPart->getClusters();
			int ntrk = tckvec.size();
			int nclu = cluvec.size();
			ntrks_Pandora += ntrk;
			nclus_Pandora += nclu;
			//
			Double_t d0=0.,z0=0.,deltad0=0.,deltaz0=0.,nsigd0=0.,nsigz0=0.;
			if (ntrk > 0) {
				d0 = tckvec[0]->getD0();
				z0 = tckvec[0]->getZ0();
				deltad0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
				deltaz0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
				nsigd0 = d0/deltad0;
				nsigz0 = z0/deltaz0;
				if( m_savehis ){
					for(int it=0; it<ntrk; it++){
						Double_t d0 = tckvec[it]->getD0();
						Double_t z0 = tckvec[it]->getZ0();
						Double_t ph = tckvec[it]->getPhi();
						//
						pd_phi.push_back(ph);
						pd_Vr .push_back(d0);
						pd_Vz .push_back(z0);
						//
						hPando[5]->Fill(ph,1.0);
						hPando[6]->Fill(d0,1.0);
						hPando[7]->Fill(z0,1.0);
					}
				}
			}
			//
			Double_t energy = recPart->getEnergy();
			Double_t charge = recPart->getCharge();
			TVector3 momentum = TVector3(recPart->getMomentum());
			Double_t cosTheta = momentum.CosTheta();
			TLorentzVector p4 = TLorentzVector(momentum,energy);
			pVis_Pandora += p4;
			if(fabs(charge)<0.01) {
				pd_pvis_neu += p4;
				neutral_pd_p4list.push_back(p4); 
			}else{
				pd_pvis_trk += p4;
				track_pd_p4list.push_back(p4); 
			}	
			if ( m_savehis ){
				hPando[10]->Fill(momentum.X(),1.0);
				hPando[11]->Fill(momentum.Y(),1.0);
				hPando[12]->Fill(momentum.Z(),1.0);
				hPando[13]->Fill(energy      ,1.0);
				hPando[14]->Fill(cosTheta    ,1.0);
				//		
				pd_costheta.push_back(cosTheta);
				pd_energy  .push_back(energy  );
			}
			if (energy>m_EnergyCut && energy<m_ECM*1.0) {
				Double_t ecalEnergy = 0;
				Double_t hcalEnergy = 0;
				Double_t yokeEnergy = 0;
				Double_t totalCalEnergy = 0;
				int nHits = 0;
				vector<lcio::Cluster*> clusters = recPart->getClusters();
				//
				for (vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster){
					ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
					hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
					yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
					ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
					hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
					CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
					nHits += calHits.size();
				}
				//
				totalCalEnergy = ecalEnergy + hcalEnergy;
				if ( m_savehis ){ 
					hPando[16]->Fill(nHits    ,1.0);
					hPando[17]->Fill(totalCalEnergy,1.0);
					pd_nHits .push_back(nHits);
					pd_calEn .push_back(totalCalEnergy);
					pd_ecalE .push_back(ecalEnergy);
					pd_hcalE .push_back(hcalEnergy);
				}

				bool woFSR = kTRUE;
				Double_t coneEnergy0[3] = {0.,0.,0.};
				Double_t pFSR[4] = {0.,0.,0.,0.};
				getConeEnergy(recPart,colPFOPandora,fCosConeCut,woFSR,coneEnergy0,pFSR);
				Double_t coneEnergy = coneEnergy0[0];
				if ( m_savehis ) hPando[18]->Fill(coneEnergy,1.0);
				Double_t coneEN     = coneEnergy0[1];
				if ( m_savehis ) hPando[19]->Fill(coneEN    ,1.0);
			}
		}
		//
		if ( m_savehis ){
			hPando[ 0]->Fill(pVis_Pandora.E (),1.0);
			hPando[ 1]->Fill(pVis_Pandora.Px(),1.0);
			hPando[ 2]->Fill(pVis_Pandora.Py(),1.0);
			hPando[ 3]->Fill(pVis_Pandora.Pz(),1.0);
			hPando[ 4]->Fill(pVis_Pandora.Pt(),1.0);
			hPando[ 8]->Fill(    ntrks_Pandora,1.0);
			hPando[ 9]->Fill(    nPFOs_Pandora,1.0);
			m_ntVal->fillDouble("pd_ntrks",       ntrks_Pandora);
			m_ntVal->fillDouble("pd_nclus",       nclus_Pandora);
			m_ntVal->fillDouble("pd_nPFOs",       nPFOs_Pandora);
			m_ntVal->fillArray ("pd_phi"         , "PD_trk", pd_phi         , pd_phi.size());
			m_ntVal->fillArray ("pd_Vz"          , "PD_trk", pd_Vz          , pd_Vz .size());
			m_ntVal->fillArray ("pd_Vr"          , "PD_trk", pd_Vr          , pd_Vr .size());
			m_ntVal->fillArray ("pd_costheta"    , "PD_pfo", pd_costheta    , pd_costheta.size());
			m_ntVal->fillArray ("pd_energy"      , "PD_pfo", pd_energy      , pd_energy.size());
			m_ntVal->fillArray ("pd_nHits"       , "PD_pfo", pd_nHits       , pd_nHits.size());
			m_ntVal->fillArray ("pd_calEn"       , "PD_pfo", pd_calEn       , pd_calEn.size());
			m_ntVal->fillArray ("pd_ecalE"       , "PD_pfo", pd_ecalE       , pd_ecalE.size());
			m_ntVal->fillArray ("pd_hcalE"       , "PD_pfo", pd_hcalE       , pd_hcalE.size());
		}
	}
	//
	if( colIsoLeps ){
		int nLeps = colIsoLeps->getNumberOfElements();
		if(m_debug>5)	printf("Lepton candidates  %8d\n", nLeps);
		for (int i=0;i<nLeps;i++) {
			ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colIsoLeps->getElementAt(i));
			//
			Double_t energy = recPart->getEnergy();
			Double_t charge = recPart->getCharge();
			if(m_debug>5)	printf("Lepton ID is %8d\n", recPart->getType());
			if ( fabs(charge)<0.01 ) continue;
			if ( energy>m_EnergyCut  && energy<m_ECM*1.0 ) {
				if( m_full>0 ){	
					int islepton = 0; 
					if ( m_useCalo == 0) islepton = abs(recPart->getType());
					else                 islepton = IsLepton(recPart);

					if ( fabs(charge)>0.01 ){
						if (  islepton==11 ) 
						{
							if ( charge>0 ) { epList.push_back(new FSParticle(recPart, colMCTL, partonList, "e+") ); }
							else            { emList.push_back(new FSParticle(recPart, colMCTL, partonList, "e-") ); }
						}
						if (  islepton==13 ) 
						{
							if ( charge>0 ) { mupList.push_back(new FSParticle(recPart, colMCTL, partonList, "mu+") ); }
							else            { mumList.push_back(new FSParticle(recPart, colMCTL, partonList, "mu-") ); }
						}
					}
				}else{

					if ( fabs(charge)>0.01 && energy > m_LepEnergyThreshold ){
						int pdgid=0;
						if( m_matchmc>0 ) {
							LCRelationNavigator *navMCTL   = new LCRelationNavigator(colMCTL);
							LCObjectVec vecMCTL            = navMCTL->getRelatedToObjects(recPart);
							if (vecMCTL.size() > 0) {
								for(unsigned int k=0; k< vecMCTL.size(); k++){
									MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL[k]);
									pdgid=abs(mcp->getPDG());
								}
							}
							delete navMCTL;
						}
						if(m_debug>3)	printf("recon object's mctruth is %8d\n", pdgid);
						if( pdgid==11){ 
							if ( charge>0 ) { epList. push_back(new FSParticle(recPart, colMCTL, partonList, "e+" ) ); }
							else            { emList. push_back(new FSParticle(recPart, colMCTL, partonList, "e-" ) ); }
						}
						if( pdgid==13){ 
							if ( charge>0 ) { mupList.push_back(new FSParticle(recPart, colMCTL, partonList, "mu+") ); }
							else            { mumList.push_back(new FSParticle(recPart, colMCTL, partonList, "mu-") ); }
						}
					}
				}
			}
		}

	}
	//
	if( colPFO ){
		_Pmax = 0 ;
		_Emax = 0 ;
		nPFOs = colPFO->getNumberOfElements();
		raw_PFOs.clear(); 
		raw_trackvec.clear(); 
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
			if(m_debug>4){
				for(int it=0; it<ntrk; it++) 
					printf("track  id of %3d %3d is %3d\n", i, it, tckvec[it]->id());
				for(int it=0; it<nclu; it++) 
					printf("clusterid of %3d %3d is %3d\n", i, it, cluvec[it]->id());
			}
			//
			Double_t d0=0.,z0=0.,deltad0=0.,deltaz0=0.,nsigd0=0.,nsigz0=0.;
			if (ntrk > 0 && m_full>0) {
				d0 = tckvec[0]->getD0();
				z0 = tckvec[0]->getZ0();
				deltad0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
				deltaz0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
				nsigd0 = d0/deltad0;
				nsigz0 = z0/deltaz0;
				//
				for(int it=0; it<ntrk; it++){
					raw_trackvec.push_back(tckvec[it]);
				}	
				//
				if( m_savehis ){
					for(int it=0; it<ntrk; it++){
						Double_t d0 = tckvec[it]->getD0();
						Double_t z0 = tckvec[it]->getZ0();
						Double_t ph = tckvec[it]->getPhi();
						//
						ab_phi.push_back(ph);
						ab_Vr .push_back(d0);
						ab_Vz .push_back(z0);
						//
						hValid[5]->Fill(ph,1.0);
						hValid[6]->Fill(d0,1.0);
						hValid[7]->Fill(z0,1.0);
					}
				}
			}
			//
			Double_t charge = recPart->getCharge();
			TVector3 momentum = TVector3(recPart->getMomentum());
			Double_t cosTheta = momentum.CosTheta();
			Double_t momentumMagnitude = momentum.Mag();
			ntrks += ntrk;
			if(fabs(charge)>0.1)nclus += nclu;
			if( energy   > _Emax ) _Emax = energy  ;
			if( momentumMagnitude  > _Pmax ) _Pmax = momentumMagnitude;
			TLorentzVector p4 = TLorentzVector(momentum,energy);
			if(fabs(charge)<0.01) {
				ab_pvis_neu += p4;
				neutral_ab_p4list.push_back(TLorentzVector(momentum,energy)); 
			}else{
				ab_pvis_trk += p4;
				track_ab_p4list.push_back(TLorentzVector(momentum,energy)); 
			}	
			if ( m_savehis ){
				hValid[10]->Fill(momentum.X(),1.0);
				hValid[11]->Fill(momentum.Y(),1.0);
				hValid[12]->Fill(momentum.Z(),1.0);
				hValid[13]->Fill(energy      ,1.0);
				hValid[14]->Fill(cosTheta    ,1.0);
				//
				ab_costheta.push_back(cosTheta);
				ab_energy  .push_back(energy  );
			}
			if ( energy>m_EnergyCut  && energy<m_ECM*1.0) {
				pVis += p4;
				if( m_full>0 ){	
					Double_t ecalEnergy = 0;
					Double_t hcalEnergy = 0;
					Double_t yokeEnergy = 0;
					Double_t totalCalEnergy = 0;
					int nHits = 0;
					vector<lcio::Cluster*> clusters = recPart->getClusters();
					for (vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster){
						ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
						hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
						yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
						ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
						hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
						CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
						nHits += calHits.size();
					}
					totalCalEnergy = ecalEnergy + hcalEnergy;

					if ( m_savehis ){ 
						hValid[16]->Fill(nHits    ,1.0);
						hValid[17]->Fill(totalCalEnergy,1.0);
						ab_nHits .push_back(nHits);
						ab_calEn .push_back(totalCalEnergy);
						ab_ecalE .push_back(ecalEnergy);
						ab_hcalE .push_back(hcalEnergy);
					}
					bool woFSR = kTRUE;
					Double_t coneEnergy0[3] = {0.,0.,0.};
					Double_t pFSR[4] = {0.,0.,0.,0.};
					getConeEnergy(recPart,colPFO,fCosConeCut,woFSR,coneEnergy0,pFSR);
					Double_t coneEnergy = coneEnergy0[0];
					if ( m_savehis ) hValid[18]->Fill(coneEnergy,1.0);
					Double_t coneEN     = coneEnergy0[1];
					Double_t coneEC     = coneEnergy0[2];
					if ( m_savehis ) hValid[19]->Fill(coneEN+coneEC,1.0);
					// select the leptons and photons 
					if(m_debug>4){
						printf("total cal energy = %10.4f, p= %10.4f, ecal energy = %10.4f, yoke energy = %10.4f\n",
								totalCalEnergy, momentumMagnitude, ecalEnergy, yokeEnergy );
					}
					if ( fabs(charge)>0.01 ){
					}else{
						//if( IsPhoton(recPart,colPFO) ) photonList.push_back(new FSParticle(recPart, colMCTL, "gamma") ); 
						//if( energy > 2 &&  ecalEnergy/totalCalEnergy >0.9) photonList.push_back(new FSParticle(recPart, colMCTL, partonList, "gamma") ); 
						if( energy > 2 ) photonList.push_back(new FSParticle(recPart, colMCTL, partonList, "gamma") );
					}
				}else{
					if ( fabs(charge)>0.01 ){
					}else{
						int pdgid=0;
						if( m_matchmc>0 ) {
							LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
							LCObjectVec vecMCTL          = navMCTL->getRelatedToObjects(recPart);
							if (vecMCTL.size() > 0) {
								for(unsigned int k=0; k< vecMCTL.size(); k++){
									MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL[k]);
									pdgid=abs(mcp->getPDG());
								}
								//if( m_pdgid != pdgid ) printf("recon object is %8d, mctruth is %8d\n", pdgid, m_pdgid);
							}
							delete navMCTL;
						}
						if( pdgid==22 && energy > 2 ) photonList.push_back(new FSParticle(recPart, colMCTL, partonList, "gamma") );
					}
				}
			}
		}
		nPFOs = raw_PFOs.size(); 
		if ( m_savehis ){
			hValid[ 0]->Fill(pVis.E (),1.0);
			hValid[ 1]->Fill(pVis.Px(),1.0);
			hValid[ 2]->Fill(pVis.Py(),1.0);
			hValid[ 3]->Fill(pVis.Pz(),1.0);
			hValid[ 4]->Fill(pVis.Pt(),1.0);
			hValid[ 8]->Fill(    ntrks,1.0);
			hValid[ 9]->Fill(    nPFOs,1.0);
			//
			m_ntVal->fill4Momentum("pd_trk", "pd_trk",   track_pd_p4list,   track_pd_p4list.size());
			m_ntVal->fill4Momentum("pd_neu", "pd_neu", neutral_pd_p4list, neutral_pd_p4list.size());
			m_ntVal->fill4Momentum("ab_trk", "ab_trk",   track_ab_p4list,   track_ab_p4list.size());
			m_ntVal->fill4Momentum("ab_neu", "ab_neu", neutral_ab_p4list, neutral_ab_p4list.size());

			m_ntVal->fill4Momentum("ab_tot_neu", ab_pvis_neu ); 
			m_ntVal->fill4Momentum("ab_tot_trk", ab_pvis_trk ); 
			m_ntVal->fill4Momentum("pd_tot_neu", pd_pvis_neu ); 
			m_ntVal->fill4Momentum("pd_tot_trk", pd_pvis_trk ); 

			m_ntVal->fillDouble("ab_ntrks",       ntrks);
			m_ntVal->fillDouble("ab_nclus",       nclus);
			m_ntVal->fillDouble("ab_nPFOs",       nPFOs);
			m_ntVal->fillArray ("ab_phi"         , "ab_trk", ab_phi         , ab_phi.size());
			m_ntVal->fillArray ("ab_Vz"          , "ab_trk", ab_Vz          , ab_Vz .size());
			m_ntVal->fillArray ("ab_Vr"          , "ab_trk", ab_Vr          , ab_Vr .size());
			m_ntVal->fillArray ("ab_costheta"    , "ab_pfo", ab_costheta    , ab_costheta.size());
			m_ntVal->fillArray ("ab_energy"      , "ab_pfo", ab_energy      , ab_energy.size());
			m_ntVal->fillArray ("ab_nHits"       , "ab_pfo", ab_nHits       , ab_nHits.size());
			m_ntVal->fillArray ("ab_calEn"       , "ab_pfo", ab_calEn       , ab_calEn.size());
			m_ntVal->fillArray ("ab_ecalE"       , "ab_pfo", ab_ecalE       , ab_ecalE.size());
			m_ntVal->fillArray ("ab_hcalE"       , "ab_pfo", ab_hcalE       , ab_hcalE.size());
		}
	}
	if ( m_savehis ) {
		//
		m_ntVal->fillDouble("ab_AllMass", pVis.M() );
		m_ntVal->fillDouble("pd_AllMass", pVis_Pandora.M() );
		m_ntVal->fillDouble("HiggsDecay", _nhfs );
		m_ntVal->write();
		//
		hValid[20]->Fill(pVis.M(), 1.0);
		if( _nhbb == 2 )                                          hValid[21]->Fill(pVis.M(), 1.0);
		if( _nhcc == 2 )                                          hValid[22]->Fill(pVis.M(), 1.0);
		if( _nhgg == 2 || _nhuu == 2 || _nhdd == 2 || _nhss == 2) hValid[23]->Fill(pVis.M(), 1.0);
		if( _nhpp == 2 )                                          hValid[24]->Fill(pVis.M(), 1.0);
		if( _nhee == 2 || _nhmm == 2 )                            hValid[25]->Fill(pVis.M(), 1.0);
		if( _nhzz == 2 )                                          hValid[26]->Fill(pVis.M(), 1.0);
		if( _nhww == 2 )                                          hValid[27]->Fill(pVis.M(), 1.0);
		if( _nhtt == 2 )                                          hValid[28]->Fill(pVis.M(), 1.0);
		//
		hPando[20]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhbb == 2 )                                          hPando[21]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhcc == 2 )                                          hPando[22]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhgg == 2 || _nhuu == 2 || _nhdd == 2 || _nhss == 2) hPando[23]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhpp == 2 )                                          hPando[24]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhee == 2 || _nhmm == 2 )                            hPando[25]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhzz == 2 )                                          hPando[26]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhww == 2 )                                          hPando[27]->Fill(pVis_Pandora.M(), 1.0);
		if( _nhtt == 2 )                                          hPando[28]->Fill(pVis_Pandora.M(), 1.0);
	}
	//
	if( colTaus ){
		numberTaus = colTaus ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of tau       = %3d\n",numberTaus ) ;
		}
		if ( numberTaus > 0) {
			for (int j=0 ; j < numberTaus ; j++ ) {
				ReconstructedParticle *pTau = dynamic_cast<ReconstructedParticle *>( colTaus->getElementAt (j) );
				Double_t charge  = pTau->getCharge();
				//
				if( charge> 0.1)
					taupList.push_back(new FSParticle(pTau, colMCTL, partonList,  "tau+" ));
				if( charge<-0.1)
					taumList.push_back(new FSParticle(pTau, colMCTL, partonList,  "tau-" ));
				//
				if(m_debug>2){
					ReconstructedParticleVec vPart = pTau->getParticles();
					int nPart= vPart.size();
					printf("tau%2d  %3d has %3d objects\n", (int)charge, j, nPart);
					if(nPart>0){
						for (int k=0 ; k < nPart ; k++ ) {
							TrackVec   tckvec = vPart[k]->getTracks();
							ClusterVec cluvec = vPart[k]->getClusters();
							int ntrk = tckvec.size();
							int nclu = cluvec.size();
							if(m_debug>3){
								printf("obj %3d of tau  %3d has %3d tacks and %3d clusters\n", k, j, ntrk, nclu);
								if( m_full>0 ){	
									for(int it=0; it<ntrk; it++) 
										printf("track  id of %3d %3d is %3d\n", j, it, tckvec[it]->id());
									for(int it=0; it<nclu; it++) 
										printf("clusterid of %3d %3d is %3d\n", j, it, cluvec[it]->id());
								}
							}
						}
					}
				}
			}
		}
	}
	//
	_y12 =  -1.0; 
	_y23 =  -1.0;
	_y34 =  -1.0;
	_y45 =  -1.0;
	_y56 =  -1.0;
	_y67 =  -1.0;
	if( col2Jets ){
		numberJets = col2Jets ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of jets in Col= %3d\n",numberJets ) ;
		}
		if ( numberJets > 0) {
			double btag=-1, ctag=-1, bctag=-1, cat=-1, flavor=0;
			PIDHandler pidh(col2Jets);
			int algo   = 0;
			int ibtag  =-1; 
			int ictag  =-1; 
			int ibctag =-1;
			int icat   =-1;
			if (m_TagFlavor>0 ) {
				algo = pidh.getAlgorithmID( "lcfiplus" );
				//algo = pidh.getAlgorithmID( "LCFIFlavourTag" );
				if( m_debug>2 ) printf("algo of lcfipus is %3d\n",algo);
				ibtag     = pidh.getParameterIndex (algo,  "BTag");	
				ictag     = pidh.getParameterIndex (algo,  "CTag");	
				ibctag    = pidh.getParameterIndex (algo, "BCTag");
				icat      = pidh.getParameterIndex (algo, "Category");
			}
			//	
			for (int j=0 ; j < numberJets ; j++ ) {
				ReconstructedParticle *pJet = dynamic_cast<ReconstructedParticle *>( col2Jets->getElementAt (j) );
				TVector3 momentum = TVector3(pJet->getMomentum());
				//
				if(m_debug>2){
					ReconstructedParticleVec vPart = pJet->getParticles();
					int nPart= vPart.size();
					int ntrkj=0, ncluj=0;
					if(nPart>0){
						for (int k=0 ; k < nPart ; k++ ) {
							TrackVec   tckvec = vPart[k]->getTracks();
							ClusterVec cluvec = vPart[k]->getClusters();
							int ntrk = tckvec.size();
							int nclu = cluvec.size();
							ntrkj+=ntrk;
							ncluj+=nclu;
							if(m_debug>4){
								printf("obj %3d of jet  %3d has %3d tacks and %3d clusters\n", k, j, ntrk, nclu);
								if( m_full>0 ){	
									for(int it=0; it<ntrk; it++) 
										printf("track  id of %3d %3d is %3d\n", j, it, tckvec[it]->id());
									for(int it=0; it<nclu; it++) 
										printf("clusterid of %3d %3d is %3d\n", j, it, cluvec[it]->id());
								}
							}
						}
					}
					printf("jet  %3d has %3d(%3d+%3d) objects\n", j, nPart, ntrkj, ncluj);
				}
				//
				if ( m_TagFlavor>0 ) {
					try{
						if( ibtag>=0 && ictag>=0 ){
							const ParticleID &pid = pidh.getParticleID(pJet, algo);
							btag  = pid.getParameters()[ibtag ];
							ctag  = pid.getParameters()[ictag ];
							bctag = pid.getParameters()[ibctag];
							cat   = pid.getParameters()[icat  ];
						}
						if(m_debug>4){
							printf("ibtag, ictag, ibctag = %4d, %4d, %4d, %4d\n",ibtag,ictag,ibctag,icat ) ;
							printf(" btag,  ctag,  bctag = %4.2f, %4.2f, %4.2f, %4.2f\n",btag,ctag,bctag,cat ) ;
						}
					}catch(...){
						btag=-1, ctag=-1, bctag=-1, flavor=0;
					}
				}
				jetList .push_back(new FSParticle(pJet, colMCTL, partonList,  "jet", btag, ctag, bctag, cat, flavor));
			}
			//
			if ( m_full>0 || m_TagFlavor>0 ) {
				ReconstructedParticle *Jet0 = dynamic_cast<ReconstructedParticle *>(col2Jets->getElementAt(0));
				Int_t algo_y = pidh.getAlgorithmID("yth");
				const ParticleID & ythID = pidh.getParticleID(Jet0, algo_y);
				FloatVec params_y = ythID.getParameters();
				_y12 = params_y[pidh.getParameterIndex(algo_y, "y12")];
				_y23 = params_y[pidh.getParameterIndex(algo_y, "y23")];
				_y34 = params_y[pidh.getParameterIndex(algo_y, "y34")];
				_y45 = params_y[pidh.getParameterIndex(algo_y, "y45")];
				_y56 = params_y[pidh.getParameterIndex(algo_y, "y56")];
				_y67 = params_y[pidh.getParameterIndex(algo_y, "y67")];
			}else{
				_y12 =  col2Jets->parameters().getFloatVal( "y12" );
				_y23 =  col2Jets->parameters().getFloatVal( "y23" );
				_y34 =  col2Jets->parameters().getFloatVal( "y34" );
				_y45 =  col2Jets->parameters().getFloatVal( "y45" );
				_y56 =  col2Jets->parameters().getFloatVal( "y56" );
				_y67 =  col2Jets->parameters().getFloatVal( "y67" );
			}
			_y12 =  _y12>1.01 ?  _y12-100.0 : _y12;  
			_y23 =  _y23>1.01 ?  _y23-100.0 : _y23;  
			_y34 =  _y34>1.01 ?  _y34-100.0 : _y34;  
			_y45 =  _y45>1.01 ?  _y45-100.0 : _y45;  
			_y56 =  _y56>1.01 ?  _y56-100.0 : _y56;  
			_y67 =  _y67>1.01 ?  _y67-100.0 : _y67;  
		}
	}
	if( col4Jets ){
		numberJets = col4Jets ->getNumberOfElements();
		if(m_debug>1){
			printf("No. of jets in 4 jet Col= %3d\n",numberJets ) ;
		}
		if ( numberJets > 0) {
			double btag=-1, ctag=-1, bctag=-1, cat=-1, flavor=0;
			PIDHandler pidh(col4Jets);
			int algo   = 0;
			int ibtag  =-1; 
			int ictag  =-1; 
			int ibctag =-1;
			int icat   =-1;
			if (m_TagFlavor>0 ) {
				algo = pidh.getAlgorithmID( "lcfiplus" );
				//algo = pidh.getAlgorithmID( "LCFIFlavourTag" );
				if( m_debug>2 ) printf("algo of lcfipus is %3d\n",algo);
				ibtag     = pidh.getParameterIndex (algo,  "BTag");	
				ictag     = pidh.getParameterIndex (algo,  "CTag");	
				ibctag    = pidh.getParameterIndex (algo, "BCTag");
				icat      = pidh.getParameterIndex (algo, "Category");
			}
			//	
			for (int j=0 ; j < numberJets ; j++ ) {
				ReconstructedParticle *pJet = dynamic_cast<ReconstructedParticle *>( col4Jets->getElementAt (j) );
				TVector3 momentum = TVector3(pJet->getMomentum());
				//
				if ( m_TagFlavor>0 ) {
					try{
						if( ibtag>=0 && ictag>=0 ){
							const ParticleID &pid = pidh.getParticleID(pJet, algo);
							btag  = pid.getParameters()[ibtag ];
							ctag  = pid.getParameters()[ictag ];
							bctag = pid.getParameters()[ibctag];
							cat   = pid.getParameters()[icat  ];
						}
						if(m_debug>4){
							printf("ibtag, ictag, ibctag = %4d, %4d, %4d, %4d\n",ibtag,ictag,ibctag,icat ) ;
							printf(" btag,  ctag,  bctag = %4.2f, %4.2f, %4.2f, %4.2f\n",btag,ctag,bctag,cat ) ;
						}
					}catch(...){
						btag=-1, ctag=-1, bctag=-1, flavor=0;
					}
				}
				jet4List .push_back(new FSParticle(pJet, colMCTL, partonList,  "jet", btag, ctag, bctag, cat, flavor));
			}
		}
	}
}



double  FSClasserProcessor::getMassColApp(TLorentzVector p4plus, TLorentzVector p4minus, double missingx, double missingy  ){

	double a = ( missingy*p4minus.X() - missingx*p4minus.Y() ) / ( p4plus.Y()*p4minus.Y() -  p4plus.X()*p4minus.Y() );
	double b = ( missingy* p4plus.X() - missingx* p4plus.Y() ) / ( p4plus.X()*p4minus.Y() -  p4plus.Y()*p4minus.X() );
	TLorentzVector nufortau1 = TLorentzVector( a* p4plus.X(), a* p4plus.Y(), a* p4plus.Z(), a* p4plus.Vect().Mag() );
	TLorentzVector nufortau2 = TLorentzVector( b*p4minus.X(), b*p4minus.Y(), b*p4minus.Z(), b*p4minus.Vect().Mag() );

	double mass_colapp   = ( p4plus + p4minus + nufortau1 + nufortau2 ).M();
	mass_colapp   = ( p4plus + p4minus + (m_ecms-totalP4) ).M();
	//double energy_colapp = ( p4plus + p4minus + nufortau1 + nufortau2 ).E();
	return mass_colapp;

}
//**********************************************
//**********************************************
int  FSClasserProcessor::IsLepton( ReconstructedParticle* pfo ) {
	const double _electronMinEnergyDepositByMomentum = 0.7; 
	const double _electronMaxEnergyDepositByMomentum = 1.4;
	const double _electronMinEcalToHcalFraction      = 0.9;
	const double _electronMaxEcalToHcalFraction      = 1.0;

	const double _muonMinEnergyDepositByMomentum     = 0.0;
	const double _muonMaxEnergyDepositByMomentum     = 0.3;
	const double _muonMinEcalToHcalFraction          = 0.0;
	const double _muonMaxEcalToHcalFraction          = 0.4;

	double cale[2]={0.0, 0.0};
	getCalEnergy( pfo , cale );
	double ecale  = cale[0];
	double hcale  = cale[1];
	double p      = TVector3( pfo->getMomentum() ).Mag();
	double calByP = p>0 ? (ecale + hcale)/p : 0;
	double calSum = ecale+hcale;
	double ecalFrac = calSum>0 ? ecale / calSum : 0;

	if(m_debug>3)	printf("Lepton ecal and hcal are  %6.3f and %6.3f\n", ecale, hcale);
	// electron
	if (     calByP   >= _electronMinEnergyDepositByMomentum
			&& calByP   <= _electronMaxEnergyDepositByMomentum
			&& ecalFrac >= _electronMinEcalToHcalFraction
			&& ecalFrac <= _electronMaxEcalToHcalFraction ){
		return 11;
	}

	// muon
	if (     calByP   >= _muonMinEnergyDepositByMomentum
			&& calByP   <= _muonMaxEnergyDepositByMomentum
			&& ecalFrac >= _muonMinEcalToHcalFraction
			&& ecalFrac <= _muonMaxEcalToHcalFraction ){
		return 13;
	}

	return abs(pfo->getType());

	//return 0;
}
//**********************************************
//**********************************************
double FSClasserProcessor::sphericity() {

	if( m_Sphericity<0 && raw_PFOs.size()>1){
		// assumes I'm a flat jet (no vertex structure)
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
		for (unsigned int i = 0; i<raw_PFOs.size(); i++) {
			ReconstructedParticle* v = raw_PFOs[i];
			TLorentzVector trkVec(v->getMomentum()[0],v->getMomentum()[1],v->getMomentum()[2],v->getEnergy());
			if( m_HiggsBoostVector.Mag()>0 )trkVec.Boost(m_HiggsBoostVector);
			sphMat(0,0) += trkVec.X()*trkVec.X();
			sphMat(0,1) += trkVec.X()*trkVec.Y();
			sphMat(0,2) += trkVec.X()*trkVec.Z();
			sphMat(1,0) += trkVec.Y()*trkVec.X();
			sphMat(1,1) += trkVec.Y()*trkVec.Y();
			sphMat(1,2) += trkVec.Y()*trkVec.Z();
			sphMat(2,0) += trkVec.Z()*trkVec.X();
			sphMat(2,1) += trkVec.Z()*trkVec.Y();
			sphMat(2,2) += trkVec.Z()*trkVec.Z();
		}
		double norm = sphMat(0,0)+sphMat(1,1)+sphMat(2,2);
		sphMat *= 1.5/norm;
		//
		TMatrixDSymEigen EigMat(sphMat);
		TVectorD eig   = EigMat.GetEigenValues();
		/*
		TMatrixD mat   = EigMat.GetEigenVectors();
		TMatrixD mat2  = mat.Invert()*mat;
		sphMat.Print("");
		mat .Print("");
		mat2.Print("");
		eig .Print();
		*/
		m_Sphericity =  eig(1)+eig(2);
		m_Aplanarity =  eig(2);
		m_C          =  (eig(0)*eig(1)+eig(0)*eig(2)+eig(1)*eig(2))*2;
		m_D          =  (eig(0)*eig(1)*eig(2))*2*9;
		//printf("eigen = %10.4f %10.4f %10.4f %10.4f\n",eig(0), eig(1), eig(2), eig(0)+eig(1)+eig(2));
	}
	return m_Sphericity; 
}
//**********************************************
//**********************************************
double FSClasserProcessor::Lsphericity() {

	if( m_Sphericity<0 && raw_PFOs.size()>1){
		// assumes I'm a flat jet (no vertex structure)
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

		double norm=0;
		for (unsigned int i = 0; i<raw_PFOs.size(); i++) {
			ReconstructedParticle* v = raw_PFOs[i];
			TLorentzVector trkVec(v->getMomentum()[0],v->getMomentum()[1],v->getMomentum()[2],v->getEnergy());
			if( m_HiggsBoostVector.Mag()>0 )trkVec.Boost(m_HiggsBoostVector);
			double pMag = (trkVec.Vect()).Mag();
			norm += pMag;
			sphMat(0,0) += (trkVec.X()*trkVec.X())/pMag;
			sphMat(0,1) += (trkVec.X()*trkVec.Y())/pMag;
			sphMat(0,2) += (trkVec.X()*trkVec.Z())/pMag;
			sphMat(1,0) += (trkVec.Y()*trkVec.X())/pMag;
			sphMat(1,1) += (trkVec.Y()*trkVec.Y())/pMag;
			sphMat(1,2) += (trkVec.Y()*trkVec.Z())/pMag;
			sphMat(2,0) += (trkVec.Z()*trkVec.X())/pMag;
			sphMat(2,1) += (trkVec.Z()*trkVec.Y())/pMag;
			sphMat(2,2) += (trkVec.Z()*trkVec.Z())/pMag;
		}
		sphMat *= 1.5/norm;
		TMatrixDSymEigen EigMat(sphMat);
		TVectorD eig   = EigMat.GetEigenValues();
		//
		m_Sphericity =  eig(1)+eig(2);
		m_Aplanarity =  eig(2);
		m_C          =  (eig(0)*eig(1)+eig(0)*eig(2)+eig(1)*eig(2))*2;
		m_D          =  (eig(0)*eig(1)*eig(2))*2*9;
		//printf("eigen = %10.4f %10.4f %10.4f %10.4f\n",eig(0), eig(1), eig(2), eig(0)+eig(1)+eig(2));
	}
	return m_Sphericity; 
}
//**********************************************
//**********************************************
//**********************************************
//**********************************************
vector<TVector3> p3list;
void FCN0(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	const double theta = par[0];
	const double phi   = par[1];
	TVector3 n0(sin(theta)*cos(phi),sin(theta)*sin(phi), cos(theta));
	//calculate T 
	double chisq = 0, sump=0;
	for (unsigned int i = 0; i<p3list.size(); i++) {
		TVector3 v = p3list[i];
		double  p = v.Mag();
		chisq    += fabs(v.Dot(n0));
	   sump     += p;	
	}
	//
	f = -chisq/sump;
}
//**********************************************
//**********************************************
void FCNP(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
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
	for (unsigned int i = 0; i<p3list.size(); i++) {
		TVector3 v = p3list[i];
		double  p = v.Mag();
		chisq    += fabs(v.Dot(n0));
	   sump     += p;	
	}
	//
	f = -chisq/sump + penalty*penalty*999;
}
//**********************************************
//**********************************************
double FSClasserProcessor::thrust() {
   
	if( m_Thrust<0 && raw_PFOs.size()>1){
		p3list.clear();
		TVector3 p3(0,0,0);
		for (unsigned int i = 0; i<raw_PFOs.size(); i++) {
			ReconstructedParticle* v = raw_PFOs[i];
			TLorentzVector trkVec(v->getMomentum()[0],v->getMomentum()[1],v->getMomentum()[2],v->getEnergy());
			if( m_HiggsBoostVector.Mag()>0 )trkVec.Boost(m_HiggsBoostVector);
			TVector3 v3(trkVec.Vect());
			p3list.push_back(v3);
			p3 += v3;
		}	
		//
		void FCN0(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
		void FCNP(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
		TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
		gMinuit->SetFCN(FCN0);
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
		gMinuit->mnexcm("SIMPLEX", arglist , 2, ierflg); // mimimization here ... 
		gMinuit->mnexcm("MIGRAD" , arglist , 2, ierflg); // mimimization here ... 
		//gMinuit->mnexcm("MINOS" , arglist , 2, ierflg); // Minos error 
		//Print results
		Double_t amin,edm,errdef;
		Int_t nvpar,nparx,icstat;
		gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		m_ThrEDM =  edm;
		m_Thrust =-amin;
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
		m_Theta = val[0];
		m_Phi   = val[1];
		TVector3 n0(sin(m_Theta)*cos(m_Phi), sin(m_Theta)*sin(m_Phi), cos(m_Theta));
		//printf("theta = %10.5f --> %10.5f; phi = %10.5f --> %10.5f\n", vstart[0], val[0], vstart[1], val[1]);
		delete gMinuit; gMinuit=NULL;
		//
		TVector3 pt(0,0,0);
		for (unsigned int i = 0; i<p3list.size(); i++) {
			pt += p3list[i];
		}	
		gMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 2 params
		gMinuit->SetFCN(FCNP);
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
		//gMinuit->mnexcm("MIGRAD" , arglist , 2, ierflg); // mimimization here ... 
		//gMinuit->mnexcm("MINOS" , arglist , 2, ierflg); // Minos error 
		//Print results
		gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		m_ThrEDMp =   edm;
		m_Major   = -amin;
		//
		for(Int_t p = 0; p < 2; p++)
		{ 
			gMinuit->mnpout(p, chnam, val[2+p], err[p], xlolim, xuplim, iuint);
		}
		//
		val[2] = fmod(val[2]+CLHEP::twopi+CLHEP::twopi, CLHEP::pi);
		val[3] = fmod(val[3]+CLHEP::twopi+CLHEP::twopi+CLHEP::pi, CLHEP::twopi) - CLHEP::pi;
		m_Thetap = val[2];
		m_Phip   = val[3];
		TVector3 n1(sin(m_Thetap)*cos(m_Phip),sin(m_Thetap)*sin(m_Phip), cos(m_Thetap));
		//	
		double minor = 0, sump=0;
		TVector3 nm=n1.Cross(n0);
		for (unsigned int i = 0; i<p3list.size(); i++) {
			TVector3 v = p3list[i];
			double  p = v.Mag();
			minor    += fabs(v.Dot(nm));
			sump     += p;	
		}
		m_Minor   = minor/sump;
		//
		//printf("thetap= %10.5f --> %10.5f; phip= %10.5f --> %10.5f: %11.8f\n", ng.Theta(), val[2], ng.Phi(), val[3], n0.Dot(n1));
		delete gMinuit; gMinuit=NULL;
		//
	}
	return m_Thrust; 
}
int FSClasserProcessor::NHScale( std::vector<TVector3>inputMCP, double binsizeTheta, double binsizePhi )	
{
	int       NMCP     = inputMCP.size();
	double    tmpTheta = 0; 
	double    tmpPhi   = 0; 	
	double    SL       = 10.0; //57.3; 
	TVector3  currP; 

	std::map<std::pair<int, int>, double> IndexPairToEnergy;
	IndexPairToEnergy.clear();
	std::pair<int, int> tmpIndex; 

	for(int k = 0; k < NMCP; k++)
	{
		currP           = inputMCP[k];
		tmpTheta        = currP.Theta();
		tmpPhi          = currP.Phi();
		tmpIndex.first  = int(tmpTheta*SL/binsizeTheta);
		tmpIndex.second = int(tmpPhi*SL*currP.Perp()/currP.Mag()/binsizePhi);
		if(IndexPairToEnergy.find(tmpIndex) == IndexPairToEnergy.end())
		{
			IndexPairToEnergy[tmpIndex] = currP.Mag();
		}
	}
	
	// cout<<"binsize "<<binsizeTheta<<" :NP "<<IndexPairToEnergy.size()<<endl;
	return IndexPairToEnergy.size(); 
}

double FSClasserProcessor::FDevent()
{
	std::vector<TVector3> p3list;
	double FD = 0; 
	if( _FD<0 && raw_PFOs.size()>1){
		for (unsigned int i = 0; i<raw_PFOs.size(); i++) {
			ReconstructedParticle* v = raw_PFOs[i];
			TLorentzVector trkVec(v->getMomentum()[0],v->getMomentum()[1],v->getMomentum()[2],v->getEnergy());
			if( m_HiggsBoostVector.Mag()>0 )trkVec.Boost(m_HiggsBoostVector);
			TVector3 v3(trkVec.Vect());
			p3list.push_back(v3); 
		}

		int NReSizeHit[5] = {0, 0, 0, 0, 0};
		int Scale[5]      = {2, 3, 4, 5, 6};
		int OriNHit       = NHScale(p3list, 1, 1); 

		for(int j = 0; j < 5; j++)
		{
			NReSizeHit[j] = NHScale(p3list, Scale[j], Scale[j]);
			FD += 0.2 * TMath::Log(double(OriNHit)/NReSizeHit[j])/TMath::Log(double(Scale[j]));
		}
		_FD = FD;
	}	

	return _FD;
}


