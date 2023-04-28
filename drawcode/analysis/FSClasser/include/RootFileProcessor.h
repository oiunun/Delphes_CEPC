#ifndef RootFileProcessor_h
#define RootFileProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: RootFileProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class RootFileProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new RootFileProcessor ; }


		RootFileProcessor() ;

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

		/** output root file name.
		 */
		std::string _outRootFile;

		int _nRun ;
		int _nEvt ;
} ;

#endif



