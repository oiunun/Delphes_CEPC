#include <iostream>
#include <sstream>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "RootFileProcessor.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

RootFileProcessor aRootFileProcessor ;


RootFileProcessor::RootFileProcessor() : Processor("RootFileProcessor") {
  
  // modify processor description
  _description = "RootFileProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value
  registerOptionalParameter( "OutputRootFile",
			     "Name of output root file",
			     _outRootFile,
			     std::string("output.root") );
}


TFile *output = 0;

void RootFileProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  stringstream out;
  out << _outRootFile << ends;
  output = new TFile(out.str().data(),"RECREATE");
}

void RootFileProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void RootFileProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...

  _nEvt ++ ;
}



void RootFileProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void RootFileProcessor::end(){ 
  
  output->Write();
  output->Close();
//   std::cout << "RootFileProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

