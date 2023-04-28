#ifndef _PFAFastSimu_hh_
#define _PFAFastSimu_hh_
#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <IMPL/LCFlagImpl.h>
#include "UTIL/LCRelationNavigator.h"
#include "IMPL/LCCollectionVec.h"   
#include "EVENT/ReconstructedParticle.h" 
#include <TNtuple.h>
#include <TObject.h>

#include <TTree.h>
#include <TFile.h>

#include "Utility.h"

using namespace Utility ;
class TTree;

class PFAFastSimu  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new PFAFastSimu ; }

		PFAFastSimu();

		~PFAFastSimu() {};

		void init();

		void processEvent( LCEvent * evtP );
		
		virtual void processRunHeader( LCRunHeader* run ) ;

		void end();

	protected:

		std::string _treeFileName;
		std::string _treeName;
		std::string _inputCollectionName ;
		std::string _recoParticleCollectionName ;
		std::string _mcTruthCollectionName ;

		LCFlagImpl Cluflag;

		int _nRun ;
		int _nEvt ;
		int _EnableSmearing;
		int _EnableEfficiency;
		int _EnableChargeSplitting; 
		int _EnableNeutralMerge; 
		int _EnableThetaScan; 
		int _overwrite;
		int _flavor;
		double _momentumCut;
		int _rejectNeutrino;
		TTree *_outputTree;

		int _Num; 
		unsigned int _eventNr;

		double _MCPEn;
		double _MCP[3]; 
		int _MCPID; 
		double _Eff, _SF; 
		double _FlatRnd; 
		double _costheta[1000000][5];
		/*
			double _P1[4], _P2[4], _P3[4], _P4[4];
			double _Pq[4][4], _RPq[4][4];
			double _InvM1[3], _InvM2[3], _RecoilM1[3], _RecoilM2[3];
			int _PDG[4];
			double  _Scale[4];
			double _m1, _m2; 
			int _PDG1, _PDG2; 
			int _evtType; 
			int _nHD, _HDPID; 
			int _sumPDG; 
			int _NOriginalQuark; 
			double _Polar1, _Polar2, _AngPlane;
			double _InvDJ, _RecoilDJ;		
			*/

		std::string _fileName;
		std::ostream *_output;
		EVENT::TrackVec                 _tracktrash;
		EVENT::ReconstructedParticleVec _particletrash;
		EVENT::ClusterVec               _clustertrash;	

};
#endif


