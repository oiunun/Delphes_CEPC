#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP / MARLIN_AIDA

#ifdef MARLIN_CLHEP  // only if CLHEP is available !


#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandEngine.h"


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "LGTrackSmearer.h"

namespace CLHEP{} 
using namespace CLHEP ;



namespace marlin{


	LGTrackSmearer::LGTrackSmearer(const std::vector<float>& resVec ){
                
		// copy the resolution vector parameters into a more structured vector 
		_resVec.clear();
		const unsigned int size = (int)(resVec.size() / ( (int)sizeof( TrackResolution)  / (int)sizeof(float) )) ;  
		//_resVec.resize( (size_t)size ); 

		int index = 0 ;

		for( unsigned int i=0 ; i <  size ; i++ ){

			float dPP    =  resVec[ index++ ] ;
			float dPT    =  resVec[ index++ ] ;
			float bF     =  resVec[ index++ ] ;
			float thMin  =  resVec[ index++ ] ;
			float thMax  =  resVec[ index++ ] ;

			_resVec.push_back( TrackResolution( dPP, dPT, bF, thMin, thMax ) );      
			//printf("%10f, %10f, %10f, %10f, %10f\n", dPP, dPT, bF, thMin, thMax );
		}
		//printf("%4d, %4d, %4d, %4d\n", (int)_resVec.size(), (int)resVec.size(), (int)sizeof(TrackResolution), size );
	}


	HepLorentzVector LGTrackSmearer::smearedFourVector( const HepLorentzVector& v, int pdgCode ){


		// find resolution for polar angle
		double costheta = cos(v.theta()) ;  
		double sintheta = sin(v.theta()) ;  
		double P        = v.vect().mag() ;

		//if( costheta > M_PI && theta < M_PI*2 )  theta = theta - M_PI; // need to transform to [0,pi] 


		double resolution = -1. , res1 = -1, res2 = -1, bf =3.5; 

		//printf("1====%10f, %10f, %10f, %10f, %10f %8d\n", res1, res2, bf, resolution, P, pdgCode );
		for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){

			if( costheta <= _resVec[i].ThMax  &&  costheta > _resVec[i].ThMin ) {
				res1 =  _resVec[i].DPP;
				bf   =  _resVec[i].BF; 
				res2 = std::min( 0.1,  _resVec[i].DPT / ( bf * sintheta*sintheta * P ) );
				resolution = pow(res1*res1+res2*res2, 0.5);
				//printf("====%10f, %10f, %10f, %10f, %10f %8d\n", res1, res2, bf, resolution, P, pdgCode );
				break ;
			}
		}

		//printf("2====%10f, %10f, %10f, %10f, %10f %8d\n", res1, res2, bf, resolution, P, pdgCode );
		HepLorentzVector sv( 0., 0. , 0., 0. ) ;

		if( resolution > - 1e-10  ) {

			// do the smearing ....


			double deltaP = RandGauss::shoot( 0.0 , P*P*resolution ) ;

			Hep3Vector n3v( v.vect() )  ;

			n3v.setMag( fabs(P + deltaP)  ) ;

			//       std::cout << " LGTrackSmearer::smearedFourVector P0: " 
			// 		<< P << " - P1 : " << P + deltaP 
			// 		<< " resolution: " << resolution 
			// 		<< std::endl ;


			// assume perfect electron and muon ID and
			// assign pion mass to everything else

			double mass = PION_MASS ; 

			if( std::abs( pdgCode ) == 11  )  // electron

				mass = ELECTRON_MASS ;

			else if( std::abs( pdgCode ) == 13  )  // muon

				mass = MUON_MASS ;

			sv.setVectM(  n3v  , mass  ) ;

		} 
		//printf("3====%10f, %10f, %10f, %10f, %10f %8d\n", res1, res2, bf, resolution, P, pdgCode );

		return sv ;

	}

}

#endif // MARLIN_CLHEP
