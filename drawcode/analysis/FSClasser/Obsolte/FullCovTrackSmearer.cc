#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP / MARLIN_AIDA

#ifdef MARLIN_CLHEP  // only if CLHEP is available !


#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandEngine.h"


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "FullCovTrackSmearer.h"

namespace CLHEP{} 
using namespace CLHEP ;



namespace marlin{


	FullCovTrackSmearer::FullCovTrackSmearer(const std::vector<float>& resVec ){

		_resVec.clear();
		const unsigned int size = (int)(resVec.size() / ( (int)sizeof( TrackResolution)  / (int)sizeof(double) )) ;  

		int index = 0 ;

		for( unsigned int i=0 ; i <  size ; i++ ){

			double sgD0   =  resVec[ index++ ] ;
			double sgD0PtA=  resVec[ index++ ] ;
			double sgD0PtB=  resVec[ index++ ] ;

			double sgZ0   =  resVec[ index++ ] ;     
			double sgZ0PtA=  resVec[ index++ ] ;     
			double sgZ0PtB=  resVec[ index++ ] ;     

			double sgT0   =  resVec[ index++ ] ;     

			double sgPhi  =  resVec[ index++ ] ;       
			double sgTheta=  resVec[ index++ ] ;       
			double sgPRel =  resVec[ index++ ] ;         

			double bF     =  resVec[ index++ ] ;
			double thMin  =  resVec[ index++ ] ;
			double thMax  =  resVec[ index++ ] ;

			_resVec.push_back( TrackResolution( 
						sgD0    ,
						sgD0PtA , 
						sgD0PtB ,
						sgZ0    ,
						sgZ0PtA ,
						sgZ0PtB ,
						sgT0    ,
						sgPhi   ,
						sgTheta ,
						sgPRel  ,
						bF, thMin, thMax ) );      
		}
		//printf("%4d, %4d, %4d, %4d\n", (int)_resVec.size(), (int)resVec.size(), (int)sizeof(TrackResolution), size );
	}


	HepLorentzVector FullCovTrackSmearer::smearedFourVector( const HepLorentzVector& v, int pdgCode ){


		// find resolution for polar angle
		double costheta = cos(v.theta()) ;  
		double sintheta = sin(v.theta()) ;  
		double P        = v.vect().mag() ;

		double resolution = -1. , res1 = -1, res2 = -1, bf =3.5; 

		//printf("1====%10f, %10f, %10f, %10f, %10f %8d\n", res1, res2, bf, resolution, P, pdgCode );
		for( unsigned int i=0 ; i <  _resVec.size() ; i++ ){

			if( costheta <= _resVec[i].ThMax  &&  costheta > _resVec[i].ThMin ) {
				res1 =  _resVec[i].SgPRel;
				bf   =  _resVec[i].BF; 
				res2 = std::min( 0.1,  _resVec[i].SgPRel / ( bf * sintheta*sintheta * P ) );
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

	void FullCovTrackSmearer::smearedTrack( const MCParticle * mcp, const int pdgCode, TrackImpl & ret){


		const HepLorentzVector v( mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], mcp->getEnergy());  
		const double costheta = cos(v.theta()) ;  
		const double sintheta = sin(v.theta()) ;  

		if ( costheta< _resVec[0].ThMax  && costheta> _resVec[0].ThMin ){ 
			ret.id();
			const double pt       = v.perp();
			const double p        = v.vect().mag() ;
			const double theta    = v.theta();
			const double phi      = v.phi();

			// find resolution for polar angle
			// compute momentum-dependent resolutions
			const double SigD0 = std::min(0.01, 
					pow(
						pow(_resVec[0].SgD0,2)+
						pow(_resVec[0].SgD0PtA /p/_resVec[0].BF, 2)/(sintheta*sintheta*sintheta), 
						0.5));

			const double SigZ0 = _resVec[0].SgZ0
				+ _resVec[0].SgZ0PtA * std::exp(-1.0 * std::abs(_resVec[0].SgZ0PtB) * pt);
			const double SigP = _resVec[0].SgPRel * p;
			// var(q/p) = (d(1/p)/dp)² * var(p) = (-1/p²)² * var(p)
			const double SigQOverP = SigP / (p * p);
			// shortcuts for other resolutions
			const double SigT0    = _resVec[0].SgT0;
			const double SigPhi   = _resVec[0].SgPhi;
			const double SigTheta = _resVec[0].SgTheta;
			// converstion from perigee d0,z0 to curvilinear u,v
			// d0 and u differ only by a sign
			const double SigU = SigD0;
			// project from z0 to the second axes orthogonal to the track direction
			const double SigV = SigZ0 * std::sin(theta);

			// draw random noise
			const float deltaD0    = SigD0    * RandGauss::shoot( 0.0 , 1.0 ) ;
			const float deltaZ0    = SigZ0    * RandGauss::shoot( 0.0 , 1.0 ) ;
			const float deltaT0    = SigT0    * RandGauss::shoot( 0.0 , 1.0 ) ;
			const float deltaPhi   = SigPhi   * RandGauss::shoot( 0.0 , 1.0 ) ;
			const float deltaTheta = SigTheta * RandGauss::shoot( 0.0 , 1.0 ) ;
			const float deltaP     = SigP     * RandGauss::shoot( 0.0 , 1.0 ) ;

			// smear the position
			Hep3Vector pos0(mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2]);
			Hep3Vector pos = pos0 + Hep3Vector(deltaD0 * std::sin(phi), deltaD0 * -std::cos(phi), deltaZ0);
			float  D0 = pow(pos[0]*pos[0]+pos[1]*pos[1],0.5); 
			float  Z0 = pos[2];  
			// smear the time
			// const double time = particle.time() + deltaT0;
			// smear direction angles phi,theta ensuring correct bounds
			const std::pair<double, double>  angles = ensureThetaBounds(phi + deltaPhi, theta + deltaTheta);
			// compute smeared direction vector
			const Hep3Vector dir(
					std::sin(angles.second) * std::cos(angles.first),
					std::sin(angles.second) * std::sin(angles.first),
					std::cos(angles.second));
			// compute smeared momentum vector
			const Hep3Vector mom = (p + deltaP) * dir;

			EVENT::FloatVec covMatrix;

			covMatrix.resize(15);

			for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
				covMatrix[icov] = 0;
			}

			covMatrix[0]  = ( SigD0     * SigD0        ); //sigma_d0^2
			covMatrix[2]  = ( SigPhi    * SigPhi       ); //sigma_phi0^2
			covMatrix[5]  = ( SigQOverP * SigQOverP    ); //sigma_omega^2
			covMatrix[9]  = ( SigZ0     * SigZ0        ); //sigma_z0^2
			covMatrix[14] = ( SigTheta  * SigTheta     )/pow(costheta,4); //sigma_tanl^2
			float  Omega = _resVec[0].BF * 3e-4 / ((p+deltaP)*sin(angles.second)*mcp->getCharge());
			float  TanLambda =  tan( angles.second) ;
			float vtx[3];
			vtx[0] = mcp->getVertex()[0] ;
			vtx[1] = mcp->getVertex()[1] ;
			vtx[2] = mcp->getVertex()[2] ;

			TrackStateImpl* ts = new TrackStateImpl( TrackState::AtIP,
					D0,
					angles.first,
					Omega,
					Z0,
					TanLambda,
					covMatrix,
					vtx) ;
			ret.addTrackState(ts);

			/*
			printf("3====%10f, %10f, %10f, %10f, %10f %8d\n", res1, res2, bf, resolution, P, pdgCode );
			*/
		}
	}

	void FullCovTrackSmearer::smearedCluster( const MCParticle * mcp, const int pdgCode, ClusterImpl & clu){
	}

}

#endif // MARLIN_CLHEP
