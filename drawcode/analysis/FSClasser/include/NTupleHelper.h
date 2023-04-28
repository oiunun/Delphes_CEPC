#ifndef NTUPLEHELPER_H
#define NTUPLEHELPER_H

using namespace std;

#include <string>
#include <sstream>
#include <map>
#include <assert.h>
#include "CLHEP/Vector/LorentzVector.h"
//#include <UTIL/LCFourVector.h>
#include "Utility.h"
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <EVENT/LCEvent.h>
#include <EVENT/ReconstructedParticle.h>
#include <lcio.h>
#include <MCTruthHelper.h>
#include <FSHelper.h>


using namespace lcio ;
using namespace std;
using namespace Utility;
using namespace CLHEP;
using namespace EVENT;
class FSParticle;
class FSInfo;
class MCTruthHelper;


// ***************************************************************
//
//  CLASS NTupleHelper
//
//    a utility to keep track of and fill trees
//
// ***************************************************************

class NTupleHelper{

	public:

		// constructor takes a Tuple object created by BOSS

		NTupleHelper(TTree* Tree, TLorentzVector p4=TLorentzVector(0,0,0,250));
		~NTupleHelper();


		// generic routine to fill a tree variable;
		//   if the variable doesn't yet exist, it is created

		void fillEvent  (LCEvent * evtP);
		void fillExtras (vector<FSParticle*> PFOs, string cat, const int type=11);
		void fillPFO    (ReconstructedParticle * pfo, FSParticle * fsp,int index, string tag, int m_full=1);
		void fillPFO    (           FSParticle * pfo, int index, string tag);
		void fillPFOs   (vector<ReconstructedParticle*> pfos, LCCollection *colMCTL=NULL);
		void fillJet    (           FSParticle * pfo, int index, string tag, int m_full=1, int m_ft = 1,  double kappa=0.3);
		void fillTracks (   vector<Track*> );
		void fillDouble (string name, double  value);
		void fillLong   (string name, int     value);
		void fillArray  (string name, string index_name, double* value, int size);
		void fillArray  (string name, string index_name, vector<double> value, int size);

		// fill specific information
		void fill4Momentum (string index_name, string tag, const vector<TLorentzVector>& p, const int size);
		void fill4Momentum (int index, string tag, const TLorentzVector& p);
		void fill4Momentum (string tag, const TLorentzVector& p);
		void fill4Momentum (int index, string tag);
		void fillMCTruthDT(MCTruthHelper* mc, FSInfo* fs );
		void fillMCTruth(MCTruthHelper* mc );
		double JetCharge(ReconstructedParticle* jet, const double kappa=0.3);

		// write the tree

		void write();



		
	private:
		bool m_bookingStage;
		
		TTree*    m_Tree;
		TLorentzVector  p4psi;

                
		map<string, double>    m_ntupleDoubleMap;
		map<string, double>    m_doubleMap;
		
		map<string, int>       m_ntupleIntMap;
		map<string, int>       m_IntMap;
		
		map<string, double*>   m_ntupleArrayMap;
		map<string, double*>   m_arrayMap;
		map<string, int    >   m_arraySize;
		
		
		const static int maxsize;
		bool   containsEntry(string name);
		string concatName(string prefix, string base, int index);
		string concatName(string prefix, string base);

};

#endif
