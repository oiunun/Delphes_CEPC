#ifndef PlotNtupleProcessor_h
#define PlotNtupleProcessor_h 1

#include <string>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <TTree.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

class PlotNtupleProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new PlotNtupleProcessor ; }

		PlotNtupleProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		StringVec _rootNames ;
		StringVec _treeNames ;
		StringVec _legendNames ;
		StringVec _variableNames ;
		string    _cut, _dir;
		int       _onlyMC;

} ;

#endif

