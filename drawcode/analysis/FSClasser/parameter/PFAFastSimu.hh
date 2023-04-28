#ifndef _PFAFastSimu_hh_
#define _PFAFastSimu_hh_
#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <IMPL/LCFlagImpl.h>
#include <TNtuple.h>
#include <TObject.h>

#include <TTree.h>
#include <TFile.h>

class TTree;

class PFAFastSimu  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new PFAFastSimu ; }

		PFAFastSimu();

		~PFAFastSimu() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:

		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;

		LCFlagImpl Cluflag;

		int _EnableSmearing;
		int _EnableEfficiency;
		int _EnableChargeSplitting; 
		int _EnableNeutralMerge; 
		int _EnableThetaScan; 
		int _overwrite;
		int _flavor;
		TTree *_outputTree;

		int _Num; 
		unsigned int _eventNr;

		float _MCPEn;
		float _MCP[3]; 
		int _MCPID; 
		float _Eff, _SF; 
		float _FlatRnd; 
		float _costheta[1000000][5];
		/*
			float _P1[4], _P2[4], _P3[4], _P4[4];
			float _Pq[4][4], _RPq[4][4];
			float _InvM1[3], _InvM2[3], _RecoilM1[3], _RecoilM2[3];
			int _PDG[4];
			float  _Scale[4];
			float _m1, _m2; 
			int _PDG1, _PDG2; 
			int _evtType; 
			int _nHD, _HDPID; 
			int _sumPDG; 
			int _NOriginalQuark; 
			float _Polar1, _Polar2, _AngPlane;
			float _InvDJ, _RecoilDJ;		
			*/

		std::string _fileName;
		std::ostream *_output;

};
#endif


