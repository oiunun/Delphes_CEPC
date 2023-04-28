/* LGCalibration.cc
 * $Id: LGCalibration.cc,v 1.1 2006-02-24 11:01:06 gaede Exp $ 
 */

#include "LGCalibration.h"

LGCalibration::LGCalibration(int nevt, int n1, int n2, int n3, 
	      double  en1, double  en2,double  en3, double enr ){
  obj()->setIntVal(0, nevt);
  obj()->setIntVal(1, n1);
  obj()->setIntVal(2, n2);
  obj()->setIntVal(3, n3);
    
  obj()->setDoubleVal(0, en1);
  obj()->setDoubleVal(1, en2);
  obj()->setDoubleVal(2, en3);
  obj()->setDoubleVal(3, enr);
}

int LGCalibration::getNEvt() {return obj()->getIntVal(0);}
int LGCalibration::getNhit1(){return obj()->getIntVal(1);}
int LGCalibration::getNhit2(){return obj()->getIntVal(2);}
int LGCalibration::getNhit3(){return obj()->getIntVal(3);}
double LGCalibration::getEnr1(){return obj()->getDoubleVal(0);}
double LGCalibration::getEnr2(){return obj()->getDoubleVal(1);}
double LGCalibration::getEnr3(){return obj()->getDoubleVal(2);}
double LGCalibration::getEnr4(){return obj()->getDoubleVal(3);}

