#ifndef SimpleMCEVENT_H
#define SimpleMCEVENT_H 1

#include <string>
#include <map>
#include <set>
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCIO.h"

#include "SimpleMCParticle.h"

class SimpleMCEvent  {

	public: 
		SimpleMCEvent( LCCollection* mcParticleCol ) ;
		~SimpleMCEvent() ; 

		int getRunNumber() const ;

		int getEventNumber() const ;

		double getWeight() const  ;

		void setRunNumber( int rn ) ;

		void setEventNumber( int en ) ;

		void setWeight(double w) ;
		
		void Print(int print) ;

	protected: 
		
		vector<SimpleMCParticle*>  _SimplifiedList;
		vector<SimpleMCParticle*>  _OriginnalList;
		
		bool  neglectISR;
		double  _weight ;
		double  _Ecms   ;


		int     _run    ;
		int     _event  ;

		int     _ntrk   ;
		int     _nneu   ;
		int     _nISR   ;
		int     _nFoto  ;
		int     _nH     ;
		int     _nZ     ;
		int     _nW     ;
		int     _nHww   ;
		int     _nHzz   ;
		int     _nHgg   ;
		int     _nHGG   ;

		int     _nem    ;
		int     _nep    ;
		int     _nmum   ;
		int     _nmup   ;
		int     _ntaum  ;
		int     _ntaup  ;
		int     _nnv    ;
		
		int     _nt     ;
		int     _ntbar  ;
		int     _nb     ;
		int     _nbar   ;
		int     _nc     ;
		int     _ncbar  ;
		int     _ns     ;
		int     _nsbar  ;
		int     _nu     ;
		int     _nubar  ;
		int     _nd     ;
		int     _ndbar  ;
		int     _ngluon ;

		int     _nCCbar ;
		int     _nBBbar ;
		int     _nQQbar ;
		int     _nDmes  ;
		int     _nDsmes ;
		int     _nBmes  ;
		int     _nBcmes ;
		int     _nBsmes ;

}; 


#endif /* ifndef SimpleMCEVENT__H */
