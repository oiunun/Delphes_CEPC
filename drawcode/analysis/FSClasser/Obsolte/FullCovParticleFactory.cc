#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP

#ifdef MARLIN_CLHEP  // only if CLHEP is available !

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/TrackImpl.h"

#include "FullCovParticleFactory.h"
#include <iostream>

#if ! LCIO_PATCHVERSION_GE( 1,51,3 )
#define USE_CLHEP // to activate code from LCIO header <UTIL/LCFourVector.h>
#endif

#include "UTIL/LCFourVector.h"

using namespace lcio ;


namespace marlin{

	FullCovParticleFactory::FullCovParticleFactory() :
		_momentumCut( 0.001 )  {

			_smearingVec.resize( NUMBER_OF_FASTMCPARTICLETYPES ) ;
		}


	void FullCovParticleFactory::setMomentumCut( double mCut ) {
		_momentumCut = mCut ;
	}

	void FullCovParticleFactory::setSmear( bool smear ) {
		_smear = smear ;
	}
	
	void FullCovParticleFactory::setNeutrino( bool neutrino ) {
		_rejectNeutrino = neutrino ;
	}

	void FullCovParticleFactory::registerFullCovSmearer( FullCovSmearer* sm , 
			FastMCParticleType type ) {
		_smearingVec[ type ] = sm ;

	}


	FastMCParticleType FullCovParticleFactory::getParticleType( const lcio::MCParticle* mcp ) {


		// assumes that mcp is a stable particle !

		FastMCParticleType type( UNKNOWN )  ;


		double charge =  mcp->getCharge()  ;

		if( charge > 1e-10 || charge < -1e-10  ){  

			type = CHARGED ;

		} else if(  mcp->getPDG() == 22 )  { // photon

			type = PHOTON ;

		} else if(  std::abs( mcp->getPDG() ) == 12 || std::abs( mcp->getPDG() ) == 14 ||
				std::abs( mcp->getPDG() ) == 16 || std::abs( mcp->getPDG() ) == 18 )  { // neutrinos - 18 is tau-prime

			type = NEUTRINO ;


		} else {  // treat everything else neutral hadron  

			type = NEUTRAL_HADRON ;
		}

		return type ;
	}


	ReconstructedParticle* FullCovParticleFactory::createReconstructedParticle( const MCParticle* mcp) {

		// this is where we do the fast Monte Carlo ....


#ifdef LCIO_PATCHVERSION_GE  // has been defined in 1.4.1 which has a bug fix in  LCFourVector<T>
		LCFourVector<MCParticle>  mc4V( mcp ) ;
#else
		HepLorentzVector mc4V( mcp->getMomentum()[0], mcp->getMomentum()[1],
				mcp->getMomentum()[2], mcp->getEnergy() )  ;
#endif



		FastMCParticleType type = getParticleType(mcp ) ;


		FullCovSmearer* sm = _smearingVec[ type ] ;

		int pdgid= abs(mcp->getPDG());
		if(  _rejectNeutrino && ( pdgid==12|| pdgid==14||pdgid==16) ){
			return 0 ;
		}
      
		TrackImpl  track ;
		ReconstructedParticleImpl* rec = NULL ;
		HepLorentzVector reco4v(0.,0.,0.,0.)  ;
		if( _smear && 	mc4V.vect().mag() >  _momentumCut  )  
		{	
			rec   =  new ReconstructedParticleImpl ;
			if( fabs(mcp->getCharge())>0.01 ) {
				float  p[3] ;
				sm->smearedTrack( mcp , mcp->getPDG(), track ) ;
				const TrackState *ts = track.getTrackState(TrackState::AtIP);

				float pt   =  3. * 3e-4 / std::abs( ts->getOmega()) ;
				float m = mcp->getMass();
				p[0] =  pt * std::cos( ts->getPhi() ) ;
				p[1] =  pt * std::sin( ts->getPhi() ) ;
				p[2] =  pt * ts->getTanLambda() ;
				float e = pow(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+m*m, 0.5);

				rec->setMomentum( p ) ;
				rec->setMass( m ) ;
				rec->setEnergy( e ) ;
				rec->setCharge( mcp->getCharge() ) ;
				//std::cout<<"charge pdg mass = " <<mcp->getCharge()<<" "<<mcp->getPDG()<<" "<<mcp->getMass()<<std::endl;

				float vtx[3] ;
				vtx[0] = mcp->getVertex()[0] ;
				vtx[1] = mcp->getVertex()[1] ;
				vtx[2] = mcp->getVertex()[2] ;
				//rec->setReferencePoint( vtx ) ;

				rec->setType(  mcp->getPDG() ) ;
				ReconstructedParticleImpl * dummy = new ReconstructedParticleImpl() ; 
				dummy->id();
				rec->addParticle( dummy   ) ; // dummy track to make it look like a real particle !!! memory leakage
				rec->addTrack  ( &track  ) ; 
			}
			ClusterImpl * dummy = new ClusterImpl() ; 
			rec->addCluster( dummy ) ; // dummy cluster to make it look like a real  particle !!! memory leakage
		}
		//
		return  rec ;
	}

} // namespace
#endif // MARLIN_CLHEP

