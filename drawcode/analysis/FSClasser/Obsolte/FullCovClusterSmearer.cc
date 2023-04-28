#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP / MARLIN_AIDA

#ifdef MARLIN_CLHEP  // only if CLHEP is available !


#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandEngine.h"
#include <cmath>

#include "FullCovClusterSmearer.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/MCParticleImpl.h"


namespace marlin{


	FullCovClusterSmearer::FullCovClusterSmearer(const std::vector<float>& resVec ){

		// copy the resolution vector parameters into a more structured vector 

		_resVec.resize(  resVec.size() / ( sizeof( ClusterResolution)  / sizeof(double) )  ) ;  // ==3

		int index = 0 ;

		for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){

			double A     =  resVec[ index++ ] ;
			double B     =  resVec[ index++ ] ;
			double thMin =  resVec[ index++ ] ;
			double thMax =  resVec[ index++ ] ;

			_resVec[i] = ClusterResolution( A, B , thMin, thMax ) ;      
		}
	}


	HepLorentzVector FullCovClusterSmearer::smearedFourVector( const HepLorentzVector& v, const int pdgCode ){


		// find resolution for polar angle
		double theta = cos(v.theta()) ;  

		//if( theta > M_PI && theta < M_PI*2 )  theta = theta - M_PI; // need to transform to [0, pi] 

		std::pair<double,double> resolution = std::make_pair( -1., -1. ) ; 

		//printf("1====%8f %8f\n", resolution.first, resolution.second  );
		for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){

			if( theta <= _resVec[i].ThMax  &&  theta > _resVec[i].ThMin ) {
				resolution.first  =  _resVec[i].A ;
				resolution.second =  _resVec[i].B ;
				break ;
			}
		}
		HepLorentzVector sv( 0., 0. , 0., 0. ) ;

		//printf("2====%8f %8f\n", resolution.first, resolution.second  );
		if( resolution.first > - 1e-10 ) {

			// do the smearing ....

			double E = v.e() ;

			double Eres  = std::sqrt( resolution.first * resolution.first + 
					resolution.second * resolution.second / E  )  ;

			double deltaE = RandGauss::shoot( 0.0 , E*Eres ) ;


			// assume massless clusters ...

			Hep3Vector n3v( v.vect() )  ;

			n3v.setMag( fabs(E + deltaE) ) ;

			double mass = 0. ;

			sv.setVectM(  n3v  , mass  ) ;

		} 
		//printf("3====%8f %8f %8d %8f\n", resolution.first, resolution.second, pdgCode, v.e()  );

		return sv ;

	}

	void    FullCovClusterSmearer::smearedTrack  ( const MCParticle * mcp, const int pdgCode, TrackImpl   & trk ){
	}
	void    FullCovClusterSmearer::smearedCluster( const MCParticle * mcp, const int pdgCode, ClusterImpl & clu ){
	}
}

#endif // MARLIN_CLHEP

