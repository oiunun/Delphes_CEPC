#ifndef MCTRUTHHELPER_H
#define MCTRUTHHELPER_H
#include <map>
#include <vector>
#include <string>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCFourVector.h>
#include <IMPL/MCParticleImpl.h>

using namespace std;
using namespace EVENT;


// ***************************************************************
//
//  CLASS McTruthHelper
//
//    Parse truth information using the MCParticleCol object
//      provided by BOSS
//
// ***************************************************************


class MCTruthHelper{

	public:

		// constructor
		MCTruthHelper( LCCollection * mcParticleCol);
		~MCTruthHelper();


		// return a list of MCParticle pointers sorted like Lambda, ALambda, e+,...,pi-,pi0
		vector<MCParticleImpl*> partonList         () { return m_partonList;      } 
		vector<MCParticleImpl*> AllParticleList    () { return m_allMcPList;      } 
		vector<MCParticleImpl*> finalParticleList  () { return m_finalPList;      } 
		vector<MCParticleImpl*> MCDecayParticleList() { return m_DecayList;       } 
		map<MCParticleImpl*, TLorentzVector> Rawp4List(){return m_Rawp4List ;   }
		map<MCParticleImpl*, TLorentzVector> SmearList(){return m_SmearList ;   }
		map<MCParticleImpl*, double        > PDGTGList(){return m_PDGTGList ;   }
		int                     nDecayParticles    () { return m_DecayList.size();}
		//
		vector<MCParticleImpl*>   GluonPList      ()  { return m_GluonPList      ;} 
		vector<MCParticleImpl*>   UpPList         ()  { return m_UpPList         ;} 
		vector<MCParticleImpl*>   AntiUpPList     ()  { return m_AntiUpPList     ;}
		vector<MCParticleImpl*>   DownPList       ()  { return m_DownPList       ;}
		vector<MCParticleImpl*>   AntiDownPList   ()  { return m_AntiDownPList   ;}
		vector<MCParticleImpl*>   StrangePList    ()  { return m_StrangePList    ;}
		vector<MCParticleImpl*>   AntiStrangePList()  { return m_AntiStrangePList;}
		vector<MCParticleImpl*>   CharmPList      ()  { return m_CharmPList      ;}
		vector<MCParticleImpl*>   AntiCharmPList  ()  { return m_AntiCharmPList  ;}
		vector<MCParticleImpl*>   BottomPList     ()  { return m_BottomPList     ;}
		vector<MCParticleImpl*>   AntiBottomPList ()  { return m_AntiBottomPList ;}
		vector<MCParticleImpl*>   TopPList        ()  { return m_TopPList        ;}
		vector<MCParticleImpl*>   AntiTopPList    ()  { return m_AntiTopPList    ;}
		vector<MCParticleImpl*>   ElectronPList   ()  { return m_ElectronPList   ;}
		vector<MCParticleImpl*>   PositronPList   ()  { return m_PositronPList   ;}
		vector<MCParticleImpl*>   MuonPList       ()  { return m_MuonPList       ;}
		vector<MCParticleImpl*>   AntiMuonPList   ()  { return m_AntiMuonPList   ;}
		vector<MCParticleImpl*>   TauPList        ()  { return m_TauPList        ;}
		vector<MCParticleImpl*>   AntiTauPList    ()  { return m_AntiTauPList    ;}
		//
		int getnFSRGamma   (){ return m_nFSRGamma         ;} 
		int getnPhoton     (){ return m_nPhoton           ;} 
		int getnFinalPhoton(){ return m_nFinalPhoton      ;} 
		int getnHiggs      (){ return m_nHiggs            ;} 
		int getnZ0         (){ return m_nZ0               ;} 
		int getnWp         (){ return m_nWp               ;} 
		int getnWm         (){ return m_nWm               ;} 

		int getnElectron   (){ return m_nElectron         ;} 
		int getnMuon       (){ return m_nMuon             ;} 
		int getnTau        (){ return m_nTau              ;} 
		int getnNuE        (){ return m_nNuE              ;} 
		int getnNuMu       (){ return m_nNuMu             ;} 
		int getnNuTau      (){ return m_nNuTau            ;} 

		int getnPi         (){ return m_nPi               ;} 
		int getnKaon       (){ return m_nKaon             ;} 
		int getnProton     (){ return m_nProton           ;} 
		int getnNeutron    (){ return m_nNeutron          ;}

		int getnBquark     (){ return m_nBottom           ;} 
		int getnCquark     (){ return m_nCharm            ;} 
		int getnLquark     (){ return m_nLight            ;} 
		bool FindHiggsMother(MCParticle *mcp);

		// get the total MC energy (to check that all particles are accounted for)
		double  PDGTagging(double PDGTruth);
		double  MCTotalEnergy();
		double  MCMissingEnergy();
		double  MomentumReso(MCParticleImpl * a_MCP);
		map<MCParticleImpl*, TLorentzVector> SmearedList(); 
		map<MCParticleImpl*, TLorentzVector> OriginaList(); 
                void   TopoEdit(LCCollection * mcParticleCol);
                void   initLatexCode();
                string getLatexCode ( int );
                TLorentzVector getHiggsP4( LCCollection * mcParticleCol);
                TLorentzVector getZP4( LCCollection * mcParticleCol);

		double getE_nu_H(){return _E_nu_H;}
		double getE_nu_Z(){return _E_nu_Z;}
		// count the number of FSR gammas produced in this event
		int MCFSRGamma(){  
			if (m_nFSRGamma != -1) {}					     
			else { m_nFSRGamma = nParticles(kFSRGamma); }
			return m_nFSRGamma; 
		}

		// id of particles coming from original particle (psi(2S),etc)
		int MCDecayParticle(unsigned int i){ 
			if (i >= m_DecayList.size()) return 0;
			return m_DecayList[i]->getPDG(); 
		}

		// helper function that returns the name of a particle

		string particleType(const MCParticleImpl* mcParticle);
		string particleType(int id);
		bool hasMothers(const MCParticleImpl* mcParticle, const int mid=411);

		// helper function that prints the mcParticleCol to the screen

		void printInformation(int print=0);

		int  nHiggsDaughters ( int idParticle, int idZ=0);
		int  nHiggsFinalState( );

	private:

		bool isParticle      ( const MCParticleImpl* mcParticle, int idParticle);
		int  nParticles      ( int idParticle);
		bool hasParticle     ( int idParticle);
		void nPrimaryParticle( LCCollection * mcParticleCol);

		int nVertices(  int idParent, 
				int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);
		int dVertices(  int idParent, 
				int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);
		bool hasVertex(int idParent, 
				int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);

		int nDecays(    int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);
		bool hasDecay(  int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);

		bool hasDaughters(const MCParticle* mcParticle, 
				int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);

		bool DhasDaughters(const MCParticleImpl* mcParticle, 
				int idDaughter1,
				int idDaughter2,
				int idDaughter3  = 0,
				int idDaughter4  = 0,
				int idDaughter5  = 0,
				int idDaughter6  = 0,
				int idDaughter7  = 0,
				int idDaughter8  = 0,
				int idDaughter9  = 0,
				int idDaughter10 = 0
				);

		vector<MCParticle*> getDaughters(const MCParticleImpl* mcParticle);


		template <typename T> void FreeDelAll( T & t ) {
			for(int i=0; i<t.size(); i++) delete t[i];
			//
			T tmp; 	t.swap(tmp);
		}
		
		template <typename T> void FreeAll( T & t ) {
			T tmp; 	t.swap(tmp);
		}

		static const int kHiggs	         =      25;     
		static const int kZ0   	         =      23;     
		static const int kWp   	         =      24;     
		static const int kWm   	         =     -24;     
		static const int kPsi2S	         =  100443;     
		static const int kPsi3770        =   30443;     
		static const int kGamma	         =      22;         
		static const int kFSRGamma       =     -22;        
		static const int kCluster        =      91;         
		static const int kString         =      92;         
		static const int kHc             =   10443;      
		static const int kChic0          =   10441;      
		static const int kChic1	         =   20443;      
		static const int kChic2	         =     445; 
		static const int kChic0p         =      61;      
		static const int kChic1p         =      62;      
		static const int kChic2p         =      63; 
		static const int kJpsi           =     443;        
		static const int kEtac	         =     441;        
		static const int kPhi            =     333;        
		static const int kOmega          =     223;        
		static const int kPi0 	         =     111;        
		static const int kPip 	         =     211;        
		static const int kPim 	         =    -211;       
		static const int kRho0           =     113;        
		static const int kRhop           =     213;        
		static const int kRhom           =    -213;       
		static const int kA00            =   10111;        
		static const int kA0p            =   10211;        
		static const int kA0m            =  -10211;        
		static const int kB10            =   10113;        
		static const int kB1p            =   10213;        
		static const int kB1m            =  -10213;       
		static const int kA10            =   20113;        
		static const int kA1p            =   20213;        
		static const int kA1m            =  -20213;       
		static const int kF01370         =   10221;       
		static const int kF01500         =   50221;       
		static const int kH1             =   10223;       
		static const int kH1p            =   10333;       
		static const int kA20            =     115;        
		static const int kA2p            =     215;        
		static const int kA2m            =    -215;       
		static const int kEtaprime       =     331;        
		static const int kEta 	         =     221;        
		static const int kKs	         =     310;        
		static const int kKl	         =     130;        
		static const int kKp  	         =     321;        
		static const int kKm  	         =    -321;       
		static const int kPp  	         =    2212;       
		static const int kPm  	         =   -2212;      
		static const int kN              =    2112;       
		static const int kAntiN          =   -2112;      
		static const int kEp  	         =     -11;         
		static const int kEm  	         =      11;        
		static const int kMup 	         =     -13;         
		static const int kMum 	         =      13;        
		static const int kTaup 	         =     -15;        
		static const int kTaum 	         =      15;         
		static const int kNuE            =      12;         
		static const int kNuMu           =      14;         
		static const int kNuTau          =      16;         
		static const int kAntiNuE        =     -12;        
		static const int kAntiNuMu       =     -14;        
		static const int kAntiNuTau      =     -16;        
		static const int kF0600          = 9000221;    
		static const int kK0             =     311;        
		static const int kAntiK0         =    -311;       
		static const int kKstarp         =     323;        
		static const int kKstarm         =    -323;       
		static const int kKstar0         =     313;        
		static const int kAntiKstar0     =    -313;
		static const int kK10            =   10313;        
		static const int kAntiK10        =  -10313;
		static const int kK0star0        =   10311;        
		static const int kAntiK0star0    =  -10311;        
		static const int kK0starp        =   10321;        
		static const int kK0starm        =  -10321;        
		static const int kK1p            =   10323;        
		static const int kK1m            =  -10323;        
		static const int kLambda         =    3122;
		static const int kALambda        =   -3122;
		static const int kD0             =     421;
		static const int kD0bar          =    -421;
		static const int kDp             =     411;
		static const int kDm             =    -411;
		static const int kDsp            =     431;
		static const int kDsm            =    -431;
		static const int kDstarP         =     413;
		static const int kDstarM         =    -413;
		static const int kDstar          =     423;
		static const int kAntiDstar      =    -423;
		static const int kDsstarP        =     433;
		static const int kDsstarM        =    -433;
		static const int kPsi4040        = 9000443;
		static const int kPsi4160        = 9010443;
		static const int kPsi4415        = 9020443;
		static const int kY4260          = 9030443;
		static const int kY4360          = 9040443;
		static const int kDeltapp        =    2224;
		static const int kAntiDeltapp    =   -2224;
		static const int kSigma0         =    3212;
		static const int kAntiSigma0     =   -3212;
		static const int kSigmastarm     =    3114;
		static const int kAntiSigmastarp =   -3114;
		static const int kXi0            =    3322;
		static const int kAntiXi0        =   -3322;
		static const int kXistarp        =    3314;
		static const int kAntiXistarm    =   -3314;

		map<int,              string       >    LatexCode          ;
		map<MCParticleImpl*, TLorentzVector>    m_SmearList        ;
		map<MCParticleImpl*, TLorentzVector>    m_Rawp4List        ;
		map<MCParticleImpl*, double        >    m_PDGTGList        ;
		vector<MCParticleImpl*>                 m_allMcPList       ;
		vector<MCParticleImpl*>                 m_EditedList       ;
		vector<MCParticleImpl*>                 m_McTopPList       ;
		vector<MCParticleImpl*>                 m_hDaughters       ;
		vector<MCParticleImpl*>                 m_finalPList       ;
		vector<MCParticleImpl*>                 m_partonList       ;
		vector<MCParticleImpl*>                 m_DecayList        ;
		//
		vector<MCParticleImpl*>                 m_GluonPList       ;
		vector<MCParticleImpl*>                 m_UpPList          ;
		vector<MCParticleImpl*>                 m_AntiUpPList      ;
		vector<MCParticleImpl*>                 m_DownPList        ;
		vector<MCParticleImpl*>                 m_AntiDownPList    ;
		vector<MCParticleImpl*>                 m_StrangePList     ;
		vector<MCParticleImpl*>                 m_AntiStrangePList ;
		vector<MCParticleImpl*>                 m_CharmPList       ;
		vector<MCParticleImpl*>                 m_AntiCharmPList   ;
		vector<MCParticleImpl*>                 m_BottomPList      ;
		vector<MCParticleImpl*>                 m_AntiBottomPList  ;
		vector<MCParticleImpl*>                 m_TopPList         ;
		vector<MCParticleImpl*>                 m_AntiTopPList     ;
		//
		vector<MCParticleImpl*>                 m_ElectronPList    ;
		vector<MCParticleImpl*>                 m_PositronPList    ;
		vector<MCParticleImpl*>                 m_MuonPList        ;
		vector<MCParticleImpl*>                 m_AntiMuonPList    ;
		vector<MCParticleImpl*>                 m_TauPList         ;
		vector<MCParticleImpl*>                 m_AntiTauPList     ;
		//
		int    m_nDecays       ;
		int    m_nFSRGamma     ;

		int    m_nHiggs        ;
		int    m_nZ0           ;
		int    m_nWp           ;
		int    m_nWm           ;
		
		int    m_nUp           ;
		int    m_nDown         ;
		int    m_nLight        ;
		int    m_nStrange      ;
		int    m_nCharm        ;
		int    m_nBottom       ;
		int    m_nTop          ;

		int    m_nElectron     ;
		int    m_nMuon         ;
		int    m_nTau          ;
		int    m_nNuE          ;
		int    m_nNuMu         ;
		int    m_nNuTau        ;

		int    m_nPi           ;
		int    m_nPi0          ;
		int    m_nKaon         ;
		int    m_nKs           ;
		int    m_nKl           ;
		int    m_nProton       ;
		int    m_nNeutron      ;

		int    m_nPhoton       ;
		int    m_nFinalPhoton  ;

		double m_TotalEnergy   ;
		double m_MissingEnergy ;
		
		double _E_nu_H, _E_nu_Z;
		static const double Pb[3] ;	//last two number put by hand
		static const double Pc[3] ;
		static const double Pg[3] ;
};
#endif
