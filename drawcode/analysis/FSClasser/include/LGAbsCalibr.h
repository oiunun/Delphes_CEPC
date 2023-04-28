#ifndef LGAbsCalibr_h
#define LGAbsCalibr_h 1

#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

/**            === LGAbsCalibr ==== <br>
 *  Processor makes:<br>
 *
 *      Output file for Absolute Energy Calibration <br>
 *      Create collection of energies of calorimeters <br>
 *
 *    @author V. L. Morgunov, A Zhelezov (DESY/ITEP)<br>
 *
 *
 */

namespace marlin {
  class LGAbsCalibr : public Processor {
  
  public:
  
    virtual Processor*  newProcessor() { return new LGAbsCalibr ; }
    LGAbsCalibr() ;
    //-----------------------------------------------------------------------
    virtual void init() ;
    virtual void processRunHeader( LCRunHeader* run ) ;
    virtual void processEvent( LCEvent * evt ) ; 
    virtual void check( LCEvent * evt ) ; 
    virtual void end() ;
    //-----------------------------------------------------------------------

  protected:

    int _nRun ;
    int _nEvt ;

    enum {
      ECAL1=0,
      ECAL2,
      HCAL
    };

    vector<int> _nlayer;
    vector<float> _coeff;
    vector<float> _cuts;
  } ;
} //namespace marlin
#endif



