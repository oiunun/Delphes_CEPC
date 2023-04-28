#ifndef FullCovProcessor_h
#define FullCovProcessor_h 1

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
#include "IMPL/TrackStateImpl.h"
#include "IMPL/LCFlagImpl.h"
#include <EVENT/LCParameters.h>

#include "lcio.h"
#include <string>
#include <vector>
#include "TH1D.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandEngine.h"
#include "Utility.h"
#include "cnpy.h"
#include <NTupleHelper.h>


using namespace std ;
using namespace lcio ;
using namespace EVENT ;
using namespace Utility ;

typedef std::vector<double> DoubleVec; 
namespace marlin{

	struct TrackResolution {

		TrackResolution() :
		         SgD0    (0.0) ,
		         SgD0PtA (0.0) , 
		         SgD0PtB (0.0) ,
		         SgZ0    (0.0) ,
		         SgZ0PtA (0.0) ,
		         SgZ0PtB (0.0) ,
		         SgT0    (0.0) ,
		         SgPhi   (0.0) ,
		         SgTheta (0.0) ,
		         SgPRelA (0.0) ,
		         SgPRelB (0.0) ,
			 BF      (0.0) ,
			 ThMin   (0.0) ,
			 ThMax   (0.0) {}

		TrackResolution(
		              double sgD0,
		              double sgD0PtA,
		              double sgD0PtB,
		              double sgZ0 ,     
		              double sgZ0PtA,     
		              double sgZ0PtB,     
		              double sgT0,     
		              double sgPhi,       
		              double sgTheta,       
		              double sgPRelA,         
		              double sgPRelB,         
			      double bF,
			      double thMin, 
			      double thMax) :
		         SgD0    (sgD0    ) ,
		         SgD0PtA (sgD0PtA ) , 
		         SgD0PtB (sgD0PtB ) ,
		         SgZ0    (sgZ0    ) ,
		         SgZ0PtA (sgZ0PtA ) ,
		         SgZ0PtB (sgZ0PtB ) ,
		         SgT0    (sgT0    ) ,
		         SgPhi   (sgPhi   ) ,
		         SgTheta (sgTheta ) ,
		         SgPRelA (sgPRelA ) ,
		         SgPRelB (sgPRelB ) ,
			 BF      (bF      ) ,
			 ThMin   (thMin   ) ,
			 ThMax   (thMax   ) {}


		double SgD0  ;
		double SgD0PtA  ;
		double SgD0PtB  ;
		//
		double SgZ0  ;
		double SgZ0PtA ;
		double SgZ0PtB  ;
		//
		double SgT0 ;
		// 
		double SgPhi  ;
		double SgTheta  ;
		// Relative momentum resolution.
		double SgPRelA, SgPRelB;
                //
		double BF  ;
		double ThMin ;
		double ThMax ;
	} ;

	struct ClusterResolution {
		ClusterResolution() :
			A(0.) ,
			B(0.) ,
			dThe(0.) ,
			dPhi(0.) ,
			ThMin(0.) ,
			ThMax(0.) {}

		ClusterResolution(double a, double b, double dthe, double dphi, double thMin, double thMax) :
			A(a) ,
			B(b) ,
			dThe(dthe) ,
			dPhi(dphi) ,
			ThMin(thMin) ,
			ThMax(thMax) {}

		double A ;
		double B ;
		double dThe;
		double dPhi;
		double ThMin ;
		double ThMax ;

	} ;

	typedef std::vector<TrackResolution>   ChResVec ;
	typedef std::vector<ClusterResolution> ClResVec ;

	class FullCovProcessor : public Processor {

		public:

			/** Returns a new instance of the processor.*/
			virtual Processor*  newProcessor() { return new FullCovProcessor ; }


			FullCovProcessor() ;


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
			int createReconstructedParticle( ReconstructedParticleImpl*  rec, const MCParticle * mcp, Track* trk, Cluster * cls ) ;
			int createTrack                ( TrackImpl   * trk, TrackStateImpl* ts, const MCParticle * mcp ) ;
			int createCluster              ( ClusterImpl * cls, const MCParticle * mcp ) ;
			// Wrap a periodic value back into the nominal range.
			template <typename T>
				inline T wrap_periodic(T value, T start, T range) {
					using std::floor;
					// only wrap if really necessary
					T diff = value - start;
					return ((0 <= diff) && (diff < range))
						? value
						: (value - range * floor(diff / range));
				}

			/// Calculate the equivalent angle in the [0, 2*pi) range.
			template <typename T>
				inline T radian_pos(T x) {
					return wrap_periodic<T>(x, T(0), T(2 * M_PI));
				}

			/// Calculate the equivalent angle in the [-pi, pi) range.
			template <typename T>
				inline T radian_sym(T x) {
					return wrap_periodic<T>(x, T(-M_PI), T(2 * M_PI));
				}
			// Calculates the equivalent angles phi and theta
			// in the [-pi, pi) and [0, pi] ranges by ensuring
			// the correct theta bounds
			template <typename T>
				inline std::pair<T, T> ensureThetaBounds(T phi, T theta) {
					T tmpPhi = radian_sym(phi);

					T tmpTht = std::fmod(theta, 2 * M_PI);
					if (tmpTht < -M_PI) {
						tmpTht = std::abs(tmpTht + 2 * M_PI);
					} else if (tmpTht < 0) {
						tmpTht *= -1;
						tmpPhi += M_PI;
						tmpPhi = tmpPhi > M_PI ? tmpPhi - 2 * M_PI : tmpPhi;
					}
					if (tmpTht > M_PI) {
						tmpTht = 2 * M_PI - tmpTht;
						tmpPhi += M_PI;
						tmpPhi = tmpPhi > M_PI ? (tmpPhi - 2 * M_PI) : tmpPhi;
					}

					return std::pair<T, T>(tmpPhi, tmpTht);
				}


		protected:

			ChResVec _ChResVec ;
			ClResVec _ClResVec ;
			/**  Input collection name */
			std::string _inputCollectionName ;
			std::vector<double> dat;
			std::vector<double> tag ;
			long unsigned int Ne, Nv, Pm; 

			/**  Ouput collection names */
			std::string _recoParticleCollectionName ;
			std::string _mcTruthCollectionName ;

			/** Momentum cut in GeV */
			float _momentumCut , _Scale, _Tag ;

			/** Resolutions of charged particles */
			FloatVec _initChargedRes ;

			/** Resolutions of photons */
			FloatVec _initPhotonRes ;

			/** Resolutions of photons */
			FloatVec _initNeutralHadronRes ;

			/** The particle factory */
			//IRecoParticleFactory*    _factory ;
			vector<ReconstructedParticle*>  _ptrash;
			vector<LCCollectionVec*>        _ctrash;
			EVENT::TrackVec                 _tracktrash;
			EVENT::ClusterVec               _clustertrash;	
			LCFlagImpl Cluflag;

			int _nRun ;
			int _nEvt ;
			NTupleHelper*  m_ntp;
			TH1D * h_Momentum[22];
			TH1D * h_CosTheta[22];
			TH1D * h_Mass    [13];
			TH1D * h_Vertex  [10];
			TH1D * h_Ntrks   [ 4];
			TH1D * h_VisEn   [ 3];
			TH1D * h_ResD    [ 6];
			int    m_luxury, m_saveNPZ;
			int    m_makeplots;
			int    m_smear, m_perfect, m_pid;
			int    m_rejectNeutrino;

	} ;

} // end namespace 

#endif


