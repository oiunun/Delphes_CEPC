#ifndef PhotonSelectionProcessor_h
#define PhotonSelectionProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TH1D.h"
#include <sstream>

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
 * @version $Id: PhotonSelectionProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class PhotonSelectionProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new PhotonSelectionProcessor ; }


		PhotonSelectionProcessor() ;

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

		/** Input collection name.
		 */
		std::string _colMCP ;
		std::string _colMCTL ;
		std::string _colPFOs ;
		std::string _colNewPFOs ;
		std::string _colPhotons ;
		std::string _colLeptons ;
		std::string _colZPFOs ;

		int _nRun ;
		int _nEvt ;

		TH1D *hStat;
		std::stringstream gCutName[20];
} ;

#endif



