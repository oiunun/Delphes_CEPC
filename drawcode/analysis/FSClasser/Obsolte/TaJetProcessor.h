#ifndef TaJetProcessor_h
#define TaJetProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

using namespace lcio ;
using namespace marlin ;

class TaJetProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new TaJetProcessor ; }


		TaJetProcessor() ;

		/** Called at the begin of the job before anything is read.
		 * Use to initialize the processor, e.g. book histograms.
		 */
		virtual void init() ;

		/** Called for every run.
		*/
		virtual void processRunHeader( LCRunHeader* run ) ;

		/** Called for every event - the working horse.
		*/
		virtual void processEvent( LCEvent * evt ) ; 

		virtual void check( LCEvent * evt ) ; 

		/** Called after data processing for clean up.
		*/
		virtual void end() ;


	protected:

		std::string _colMCName;
		std::string _pfoCollectionName;
		std::string _outputCollectionName;
		double _smearAngle;
		double _largestAngle;
		double _tauMass;

		int _nRun ;
		int _nEvt ;
		int _maxTracks;
} ;

#endif



