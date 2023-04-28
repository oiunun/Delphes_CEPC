#ifndef LGCalibration_h
#define LGCalibration_h 1

//C++
#include "iostream"
#include "string"

//LCIO
#include "lcio.h"
#include "UTIL/LCFixedObject.h"


#define LGCalibrationNINTVals 4  // N event and N hits in it
#define LGCalibrationNFLOATVals 0 
#define LGCalibrationNDOUBLEVals 4 // Energies

class LGCalibration : public UTIL::LCFixedObject<LGCalibrationNINTVals,
  LGCalibrationNFLOATVals,LGCalibrationNDOUBLEVals> {
  
public: 
  
    /** Convenient constructor.
     */
  LGCalibration(int nevt,  int n1,  int n2,  int n3, 
	      double en1, double  en2,double  en3, double enr );

  /** 'Copy constructor' needed to interpret LCCollection read from file/database.
   */
  LGCalibration(EVENT::LCObject* obj) : UTIL::LCFixedObject<LGCalibrationNINTVals,
							 LGCalibrationNFLOATVals,
							 LGCalibrationNDOUBLEVals>(obj) { } 
    
    /** Important for memory handling*/
  virtual ~LGCalibration(){};
  
  // the class interface:
  int getNEvt();
  int getNhit1();
  int getNhit2();
  int getNhit3();
  double getEnr1();
  double getEnr2();
  double getEnr3();
  double getEnr4();

    // -------- need to implement abstract methods from LCGenericObject
    const std::string getTypeName() const { 
	return std::string("LGCalibration");
    } 
    const std::string getDataDescription() const {
	return std::string("i:Nevent,Nhits[3],d:Energies[4]"); 
    }
    
}; // class


#endif
