#ifndef LGStdHepReader_h
#define LGStdHepReader_h 1
#include <vector>
#include "marlin/DataSourceProcessor.h"
#include "UTIL/LCStdHepRdr.h"

using namespace lcio ;

/** Reads binary StdHep files.
 *  Example processor for reading non-LCIO input files - creates events with
 *  MCParticle collections from binary StdHep files. Has to be the first active processor
 *  and requires that no LCIO input collection is used (parameter LCIOInputFiles).
 *
 *  <h4>Input - Prerequisites</h4>
 *  StdHep file.
 *
 *  <h4>Output</h4> 
 *  LCEvent with MCParticle collection.
 *
 * @param StdHepFileName   name of input file
 *
 * @author F. Gaede, DESY
 * @version $Id: LGStdHepReader.h,v 1.3 2005-10-11 12:56:28 gaede Exp $ 
 */
namespace marlin{

	class LGStdHepReader : public DataSourceProcessor {

		public:

			LGStdHepReader() ;

			virtual LGStdHepReader*  newProcessor() ;


			/** Creates events with MCParticle collections from the LGStdHep input file and
			 *  calls all active processors' processEvent() and processRunHeader Method.
			 *
			 */
			virtual void readDataSource( int numEvents ) ;

		        virtual void processMC(LCCollection* col) ;

			virtual void init() ;
			virtual void end() ;

		protected:

                        int _skip; 
			StringVec _fileNames ;
			std::vector<LCStdHepRdr*> rdr ;

	};

}
#endif
