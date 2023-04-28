#ifndef FullCovParticleFactory_h
#define FullCovParticleFactory_h 1

#include "marlin/MarlinConfig.h"

#ifdef MARLIN_CLHEP  // only if CLHEP is available !

#include "marlin/IRecoParticleFactory.h"
#include "marlin/FastMCParticleType.h"
#include "FullCovSmearer.h"


namespace marlin {

	/** Implementation of IRecoParticleFactory that implements the default behaviour 
	 *  as described in LGFastMCProcessor, i.e. have polar angle ranges with different resolutions
	 *  for charged tracks, photons and neutral hadrons.
	 *
	 *  @author F. Gaede, DESY
	 *  @version $Id: FullCovParticleFactory.h,v 1.3 2007-11-23 20:09:12 gaede Exp $ 
	 */ 

	class FullCovParticleFactory : public IRecoParticleFactory {

		public:

			FullCovParticleFactory() ;


			/** Virtual d'tor.*/
			virtual ~FullCovParticleFactory() {} 


			/** The actual factory method that creates a new ReconstructedParticle
			*/ 
			virtual lcio::ReconstructedParticle* createReconstructedParticle( const lcio::MCParticle* mcp) ;


			/** Register a particle four vector smearer for the given type.
			*/
			virtual void registerFullCovSmearer( FullCovSmearer* sm , FastMCParticleType type ) ; 


			/** Returns the type of the MCParticle. 
			*/
			virtual FastMCParticleType getParticleType( const lcio::MCParticle* mcp ) ;


			//     /** Helper function to determine the charge from the PDG (charge is missing from stdhep files) 
			//      */
			//     virtual float getCharge( int pdgCode ) ;


			/** Set the momentum cut in GeV - no particles are produced for ower momenta. Default is 0.1 eV.
			*/
			virtual void setMomentumCut( double mCut    ) ;
			virtual void setSmear      ( bool  smear    ) ;
			virtual void setNeutrino   ( bool  neutrino ) ;

		protected:

			std::vector<FullCovSmearer*> _smearingVec ;
			bool   _smear          ;
			double _momentumCut    ;
			bool   _rejectNeutrino ;

	};

}

#endif // MARLIN_CLHEP

#endif // FullCovParticleFactory_h
