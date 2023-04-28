#ifndef FSHELPER_H
#define FSHELPER_H


#include <string>
#include <map>

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
#include <TVectorD.h> 
#include <TMatrixDSym.h> 

#include "JetFitObject.h"

#include <CLHEP/Vector/LorentzVector.h>
#include "TLorentzVector.h"

#include "NTupleHelper.h"
#include "Utility.h"

using namespace Utility;
using namespace std;
using namespace lcio;
using namespace marlin;
using namespace CLHEP;

class NTupleHelper;
class FSParticle;
class FSCut;



class FSInfo{

	public:

		FSInfo(string FSName, NTupleHelper* nt, NTupleHelper* ntgen);
		~FSInfo();

		string FSName()                { return m_FSName;                               }

		bool exclusive()               { return (m_FSName.find("EXC") != string::npos); }
		bool inclusive()               { return (m_FSName.find("INC") != string::npos); }

		vector<string> particleNames() { return m_particleNames;                        }
		vector<int>    particleStatus(){ return m_particleStatus;                       }
		int nChargedParticles()        { return m_nChargedParticles;                    }
		int nFSParticles()             { return m_particleNames.size();                 }
		bool allNeutral()              { return m_nChargedParticles == 0;               }

		NTupleHelper* NT()             { return m_NT;                                   }
		NTupleHelper* NTGen()          { return m_NTGen;                                }

		int decayCode1()               { return m_decayCode1;                           }
		int decayCode2()               { return m_decayCode2;                           }
		int nMissingParticles()        { return m_nMissingParticles;                    }

		bool hasParticle(string particleName){  
			for (unsigned int i = 0; i < m_particleNames.size(); i++){
				if (m_particleNames[i] == particleName) return true; }
			return false; 
		}
	
		
		void Print(){ 
			cout << "FSClasser:  Checking the Final State " << m_FSName << endl;
			for (unsigned int i = 0; i < m_particleNames.size(); i++){
				if (m_particleStatus[i]==1)printf("FSClasser: %8s:  normal\n",m_particleNames[i].data());
				else                       printf("FSClasser: %8s: missing\n",m_particleNames[i].data());
			}
		}

		double MissingMass2(){
			if(missingMassFit()) return m_missingMassValue*m_missingMassValue; 
			else                 return 0; 
		}

		bool   constrain4Mom ()             { return m_Constrain4Mom ; }
		bool   missingMassFit()             { return m_missingMassFit; }
		double missingMassValue()           { return m_missingMassValue; }
		void   setMissingMass(double m)     { m_missingMassValue = m; }
		void   addMissingParticle(string p) {
			m_nMissingParticles++;
			m_missedParticle = p;
			m_missingMassFit = true;
			m_particleNames .push_back(p); 
			m_particleStatus.push_back(0); 
			setMissingMass(Mass(p)); 
		}


		vector< vector<FSParticle*> >& particleCombinations() { return m_particleCombinations; }


		// check to see if a certain event hypothesis satisfies the list of 
		// cuts on intermediate masses, etc.

		bool evaluateFSCuts( const vector<FSParticle*>& particleCombination,
				const TLorentzVector& pInitial, string fourVectorType );

		TVector3 HiggsBoostVector( const vector<FSParticle*>& particleCombination, string RecoilSide );

		// externally set the list of particle combinations

		void setParticleCombinations( const vector< vector<FSParticle*> >& particleCombinations)
		{ m_particleCombinations = particleCombinations; }


		// add a new FSCut (a cut on intermediate mass, for example)

		void addFSCut(FSCut* fsCut){ m_FSCuts.push_back(fsCut); }

		// static functions to unpack the final state name and parse strings

		static vector<string> getParticleNamesFromFSName   ( const string& FSName );
		static vector<int>    getDecayCodeFromParticleNames( const vector<string>& particleNames );
		static vector<int>    getDecayCodeFromFSName       ( const string& FSName );
		static vector<string> parseString                  ( const string& inputString );
		static int    getNChargedParticlesFromParticleNames( const vector<string>& particleNames, const vector<string>& particleStatus);

		// returns all combinations of indices for a submode of this final state
		// example:  final state = K+ K- pi+ pi+ pi- pi- pi0  submode = pi+ pi0
		//           returns 2,6; 3,6
		vector< vector<unsigned int> >& submodeIndices(const string& submodeName);
	 	private:
			string m_FSName;
			vector<string> m_particleNames;
			vector<int>    m_particleStatus;
			int            m_nChargedParticles;
			int            m_nMissingParticles;
			NTupleHelper*  m_NT;
			NTupleHelper*  m_NTGen;

			int m_decayCode1;
			int m_decayCode2;

			bool   m_fast;

			bool   m_Constrain4Mom ;
			bool   m_missingMassFit;
			double m_missingMassValue;
			string m_missedParticle;

			vector<FSCut*>                                m_FSCuts;
			vector<vector<FSParticle*> >                  m_particleCombinations;
			map<string, vector< vector <unsigned int> > > m_submodeIndices;


};

class FSParticle{

	public:

		FSParticle(ReconstructedParticle* pfo, LCCollection *colMCTL, 
				vector<MCParticleImpl*> partonList, string name, 
				double btag=-1, double ctag=-1, double bctag=-1, double cat =-1, double flavor=0);
		FSParticle(MCParticleImpl*            mcp, string name, TLorentzVector p4);
		FSParticle(string           name, bool missed);
		~FSParticle();
      
		// TRACK
		ReconstructedParticle *  pfo(){ return m_pfo; }
		MCParticle *             mcp(){ return m_mcp; }

		// COMMON INFORMATION
		string           name()             { return m_name;  }
		double           mass()             { return m_mass;  }
		double           energy()           { return m_Energy;}
		double           charge()           { return m_charge;}
		double           pT()               { return m_pT;}
		double           pZ()               { return m_pZ;}
		double           cosTheta()         { return m_CosTheta;}
		double           rapidity()         { return m_Rapidity;}
		double           type()             { return m_type;}
		double           btag()             { return m_btag;}
		double           cat ()             { return m_cat ;}
		double           ctag()             { return m_ctag;}
		double           bctag()            { return m_bctag;}
		double           flavor()           { return m_flavor;}
		double           leptonType()       { return m_leptonType;}
		double           sphericity(); 
		bool             missed()           { return m_missed;}
		bool             fast()             { return m_fast;  }
		void             set_miss(bool miss){ m_missed = miss;}
		JetFitObject*    jetfitobject()     { return m_JetFitObject;}
		//
		TLorentzVector rawFourMomentum()  { return m_rawFourMomentum; }
		TLorentzVector fitFourMomentum()  { return m_fitFourMomentum; }
		//
		void setRawFourMomentum (const TLorentzVector& p) {m_rawFourMomentum  = p;}
		void setFitFourMomentum (const TLorentzVector& p) {m_fitFourMomentum  = p;}


		bool duplicate(FSParticle* fsp, int full=1);
		vector<int> trackId()    { return m_trackId;    }
		vector<int> showerId()   { return m_showerId;   }
		vector<int> particleId() { return m_particleId; }

		TLorentzVector fourMomentum(ReconstructedParticle * obj );
		TLorentzVector fourMomentum(MCParticle            * obj );

		void   PrintTrackAndShowers ( );
		MCParticle * FindParton(MCParticle *mcp);
		int          FindLeptonParent(MCParticle *mcp);

	private:

		ReconstructedParticle *        m_pfo;
		MCParticle *                   m_mcp;
		JetFitObject*                  m_JetFitObject;


		string   m_name;
		int      m_type;
		int      m_pdgid;
		bool     m_missed;
		bool     m_fast;
		double   m_mass;
		double   m_recmass;
		double   m_charge;
		double   m_pT;
		double   m_pZ;
		double   m_Energy;
		double   m_Rapidity;
		double   m_CosTheta;
		double   m_btag;
		double   m_ctag;
		double   m_cat;
		double   m_bctag;
		double   m_flavor;
		double   m_leptonType;
		double   m_sphericity;

		TLorentzVector      m_rawFourMomentum;
		TLorentzVector      m_fitFourMomentum;

		vector<int>         m_trackId;
		vector<int>         m_showerId;
		vector<int>         m_particleId;
};

//********************************************************************
//********************************************************************
//********************************************************************
//
//   The FSCut class:  defines cuts on intermediate masses
//                     and recoil masses
//
//********************************************************************
//********************************************************************
//********************************************************************


class FSCut{

	public:

		FSCut(const vector<string> initialization);

		string FSName()      { return m_FSName; }
		string submodeName() { return m_submodeName; }
		string cutType()     { return m_cutType; }
		double lowCut()      { return m_lowCut;  }
		double highCut()     { return m_highCut; }

		bool Raw()     { return (m_cutType.find("Raw")     != string::npos); }
		bool Int()     { return (m_cutType.find("Int")     != string::npos); }
		bool Fit()     { return (m_cutType.find("Fit")     != string::npos); } 
		bool Recoil()  { return (m_cutType.find("Recoil")  != string::npos); }
		bool Mass()    { return (m_cutType.find("Mass")    != string::npos); }
		bool Squared() { return (m_cutType.find("Squared") != string::npos); }


	private:

		string m_FSName;
		string m_submodeName;
		string m_cutType;
		double m_lowCut;
		double m_highCut;

};
#endif
