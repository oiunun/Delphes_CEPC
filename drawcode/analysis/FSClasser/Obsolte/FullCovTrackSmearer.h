#ifndef FullCovTrackSmearer_h
#define FullCovTrackSmearer_h 1

#include "marlin/MarlinConfig.h"

#ifdef MARLIN_CLHEP  // only if CLHEP is available !

#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/MCParticleImpl.h"
#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"
#include <vector>
#include "FullCovSmearer.h"

#define ELECTRON_MASS 0.0005109989 
#define MUON_MASS     0.10565836 
#define PION_MASS     0.139570 

using namespace IMPL ;
using namespace EVENT ;

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
		         SgPRel  (0.0) ,
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
		              double sgPRel,         
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
		         SgPRel  (sgPRel  ) ,
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
		double SgPRel;
                //
		double BF  ;
		double ThMin ;
		double ThMax ;

	} ;


	class FullCovTrackSmearer : public FullCovSmearer {

		typedef std::vector<TrackResolution> ResVec ;

		public:

		FullCovTrackSmearer(const std::vector<float>& resVec ) ;


		virtual ~FullCovTrackSmearer() {} 

		virtual HepLorentzVector smearedFourVector( const HepLorentzVector& v, const int pdgCode ) ;
		virtual void     smearedTrack     ( const MCParticle    * mcp, const int pdgCode, TrackImpl   & trk ) ;
		virtual void     smearedCluster   ( const MCParticle    * mcp, const int pdgCode, ClusterImpl & clu ) ;


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

		ResVec _resVec ;

	} ;



} // end namespace 

#endif // MARLIN_CLHEP
#endif // FullCovTrackSmearer_h

