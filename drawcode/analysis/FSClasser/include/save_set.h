#ifndef _save_set_hh_
#define _save_set_hh_

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
//#include <RConfigure.h>
#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <algorithm>

#include <TLorentzVector.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>

using namespace std;
class TTree;

class save_set  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new save_set ; }

		save_set();

		~save_set() {};

		void init();

		//   void processRunHeader( LCRunHeader* run) { } 

		void processEvent( LCEvent * evtP );

		void end();
		int  FindLeptonParent(MCParticle *mcp);
		void  GetHelicityAngles( const vector<TLorentzVector> & p4List, Double_t *theta, Double_t *phi );
		void  MakePlots() ; 
		TLorentzVector HPINV(const TLorentzVector &P1, const TLorentzVector &PT);
		TLorentzVector erout4(const TLorentzVector &p1, const double phi, const double theta );

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;

		int _overwrite;
		TTree *_outputMCP; 

		std::ostream *_output;
		std::string   _fileName;
		std::string   _histFileName;

      double _Run, _Evt;
		double _CosTheta[500], _Phi[500];
		double _Mass[500];
		double _posx[500];
		double _posy[500];
		double _posz[500];
		double _posr[500];
		double _dist[500];
		double _PDGi[500];

		TH1D *h_CosTheta[6];
		TH1D *h_Phi     [6];
};



#endif


