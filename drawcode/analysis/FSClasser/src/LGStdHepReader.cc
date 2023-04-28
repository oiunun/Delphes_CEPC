
#include "LGStdHepReader.h"

#include "marlin/ProcessorMgr.h"

#include "IMPL/LCEventImpl.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/LCRunHeaderImpl.h"

#include "UTIL/LCTOOLS.h"

namespace marlin{


	LGStdHepReader aLGStdHepReader ;


	LGStdHepReader::LGStdHepReader() : DataSourceProcessor("LGStdHepReader") {

		_description = "Reads LGStdHep files as input and creates LCIO events with MCParticle collections."
			" Make sure to not specify any LCIOInputFiles in the steering in order to read LGStdHep files." ;

		registerProcessorParameter( "LGStdHepFileName" , 
				"input files"  ,
				_fileNames ,
				std::vector<std::string>() ) ;

	   registerProcessorParameter("SkipEvt","",  _skip, 0);
	}

	LGStdHepReader*  LGStdHepReader::newProcessor() { 
		return new LGStdHepReader ;
	}

	void LGStdHepReader::init() {    
		for( unsigned int ii=0; ii<_fileNames.size(); ii++)
			rdr.push_back( new  LCStdHepRdr( _fileNames[ii].c_str()) ) ;
		printParameters() ;    
	}


	void LGStdHepReader::readDataSource( int numEvents ) {

		LCCollection* col=NULL ;
		LCEventImpl*  evt=NULL ;
		static int evtNum = 0 ;
		static int runNum = 0 ;
		static unsigned int i =0; 

		while ( ((numEvents > 0 && evtNum+1 < numEvents) || numEvents==0) && i<_fileNames.size() ) {
			col = rdr[i]->readEvent();
         //if ( numEvents < _skip ) continue; 
			if ( col == NULL ){
				col = rdr[i]->readEvent();
				if ( col == NULL ){
					i++;
					runNum++; 
					_isFirstEvent = true ;	
				}
				else 
				{
					if( isFirstEvent() ) {   // create run header

						LCRunHeaderImpl* rHdr = new LCRunHeaderImpl ;

						rHdr->setDescription( " Events read from stdhep input file: " + _fileNames[i] ) ; 
						rHdr->setRunNumber( runNum ) ;

						ProcessorMgr::instance()->processRunHeader( rHdr ) ;
						_isFirstEvent = false ;	
					}

					evt = new LCEventImpl ;
					evt->setRunNumber( runNum ) ;
					evt->setEventNumber( evtNum++ ) ;

					//processMC(col);

					evt->addCollection(  col, "MCParticle"  ) ;

					ProcessorMgr::instance()->processEvent( evt ) ;

					delete evt ;
				}
			}else{

				evt = new LCEventImpl ;
				evt->setRunNumber( runNum ) ;
				evt->setEventNumber( evtNum++ ) ;

				//processMC(col);

				evt->addCollection(  col, "MCParticle"  ) ;

				ProcessorMgr::instance()->processEvent( evt ) ;

				delete evt ;
			}
		}
		//if(col!=NULL)delete col;
	}

	void LGStdHepReader::processMC(LCCollection* col){
		for(int i=0; i<col->getNumberOfElements() ; i++){

			MCParticleImpl* mcp = dynamic_cast<MCParticleImpl*> ( col->getElementAt( i ) ) ;

			if( mcp->getGeneratorStatus() == 1 ) { // stable particles only
			    double pnt[3]={100,100,100};	
             mcp->setEndpoint(pnt);
				}
			}

	}


	void LGStdHepReader::end() {

		for( unsigned int ii=0; ii<_fileNames.size(); ii++)
			if( rdr[ii] ) delete rdr[ii] ;
	}

}
