#ifndef LGFastMCProcessor_h
#define LGFastMCProcessor_h 1

#include "marlin/Processor.h"
#include "marlin/IRecoParticleFactory.h"

#include "EVENT/ReconstructedParticle.h" 
#include "UTIL/LCRelationNavigator.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ClusterImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/LCFlagImpl.h"
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>

#include "lcio.h"
#include <string>
#include <vector>
#include "TH1D.h"

#include "Utility.h"
#include <NTupleHelper.h>


using namespace std ;
using namespace lcio ;
using namespace Utility ;

namespace marlin{


	/** A simple smearing "Monte Carlo" processor.
	 *  It creates ReconstructedParticles from MCParticles according to
	 *  the resolution that is specified for the particle type, one of:
	 *  <ul>
	 *  <li>photon</li>
	 *  <li>charged: e+-,mu+-,pi+-,K+-,p+-,....</li>
	 *  <li>neutral hadron: K0L, n, Lambda0,...</li>
	 *  </ul>
	 *  The resolutions for charged particles are given as delta(1/P) for a certain polar angle 
	 *  range (mapped to [0.,pi/2.]), e.g. <br>
	 *  <b>ChargedResolution  &nbsp; .7e-5   &nbsp; 0.   &nbsp; 3.141593/2. </b><br>
	 *  sets the resolution for all charged particles to 0.7*10^-5.<br> 
	 *  Specifying different resolutions for different polar angle ranges allows to mimic degrading 
	 *  detector performance in the very forward region. <br>
	 *  The energy of the charged ReconstructedParticles is set from the "measured"
	 *  momentum using the energy momentum relation for electrons and muons and assuming the pion mass
	 *  for all other charged tracks.<br>
	 *  The resolutions for neutral particles are given as A and B for a certain polar angle range,
	 *  where dE/E = A "+" B / sqrt(E/GeV), e.g. <br>
	 *  <b>NeutralHadronResolution   &nbsp; 0.03  &nbsp; .30  &nbsp;  0.  &nbsp;  3.141593/2. </b><br>
	 *  sets the resolution for all neutral hadrons to 3% "+" 30% / sqrt( E /GeV ).
	 *  The resolution for gammas is specified in <b>PhotonResolution</b><br> 
	 *  No ReconstructedParticles are created if there is no resolution defined at that polar angle, e.g.<br>
	 *  <b>PhotonResolution  &nbsp; .7e-5   &nbsp; 0.083   &nbsp; 3.141593/2. </b><br>
	 *  effectively limits the acceptance region for photons to theta > 83mrad.<br>
	 *  
	 *  A collection of LCRelations, called "MCTruthMapping" holds the relation between the 
	 *  ReconstructedParticles and their proper MCParticles.
	 *
	 * 
	 *  <h4>Input - Prerequisites</h4>
	 *  A collection of MCParticles (the MCPArticle collection).
	 *
	 *  <h4>Output</h4> 
	 *  <ul>
	 *  <li><b>ReconstructedParticles</b>: the collection of reconstructed particles</li>
	 *  <li><b>MCTruthMapping</b>:  holds LCRelations  that map the  ReconstructedParticles to their
	 *       proper MCParticles </li>
	 *  </ul>
	 * 
	 * @param ChargedResolution    Resolution of charged particles in polar angle range:  d(1/P)  th_min  th_max
	 * @param InputCollectionName  Name of the MCParticle input collection
	 * @param MomentumCut          No reconstructed particles are produced for smaller momenta (in [GeV])
	 * @param NeutralHadronResolution Resolution dE/E=A+B/sqrt(E/GeV) of neutral hadrons in polar angle range: A  B th_min  th_max
	 * @param PhotonResolution   Resolution dE/E=A+B/sqrt(E/GeV) of photons in polar angle range: A  B th_min  th_max
	 *
	 * @param RecoParticleCollectionName    default is "ReconstructedParticles"
	 * @param MCTruthMappingCollectionName  default is "MCTruthMapping"
	 * 
	 *  @author F. Gaede, DESY
	 *  @version $Id: LGFastMCProcessor.h,v 1.4 2007-07-04 12:13:06 gaede Exp $ 
	 */

	class LGFastMCProcessor : public Processor {

		public:

			/** Returns a new instance of the processor.*/
			virtual Processor*  newProcessor() { return new LGFastMCProcessor ; }


			LGFastMCProcessor() ;


			/** Initializes ...
			*/
			virtual void init() ;

			/** Called for every run.
			*/
			virtual void processRunHeader( LCRunHeader* run ) ;

			/** Updates all registered conditions handlers and adds the data to the event.
			*/
			virtual void processEvent( LCEvent * evt ) ; 

			/** Creates some checkplots.
			*/
			virtual void check( LCEvent * evt ) ; 

			/** Called after data processing for clean up.
			*/
			virtual void end() ;

			int FindParton(MCParticle *mcp);

		protected:

			/**  Input collection name */
			std::string _inputCollectionName ;

			/**  Ouput collection names */
			std::string _recoParticleCollectionName ;
			std::string _mcTruthCollectionName ;

			/** Momentum cut in GeV */
			float _momentumCut ;

			/** Resolutions of charged particles */
			FloatVec _initChargedRes ;

			/** Resolutions of photons */
			FloatVec _initPhotonRes ;

			/** Resolutions of photons */
			FloatVec _initNeutralHadronRes ;

			/** The particle factory */
			IRecoParticleFactory*    _factory ;
			vector<ReconstructedParticle*>  _ptrash;
			vector<LCCollectionVec*>        _ctrash;
			EVENT::TrackVec                 _tracktrash;
			EVENT::ClusterVec               _clustertrash;	
			LCFlagImpl Cluflag;

			int _nRun ;
			int _nEvt ;
			NTupleHelper*  m_ntp;
			TH1D * h_Momentum[30];
			TH1D * h_CosTheta[20];
			TH1D * h_Mass    [10];
			TH1D * h_Vertex  [10];
			TH1D * h_Ntrks   [ 4];
			TH1D * h_VisEn   [ 3];
			int    m_luxury;
			int    m_makeplots;
			int    m_smear;
			int    m_rejectNeutrino;

	} ;

} // end namespace 

#endif


