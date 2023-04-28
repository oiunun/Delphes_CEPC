/**
 * @brief Marlin processor for finding isolated leptons.
 * @author Ryo Yonamine <yonamine@post.kek.jp>
 * @author Tomohiko Tanabe <tomohiko@icepp.s.u-tokyo.ac.jp>
 * @version $Id:$
 *
 * Given a list of ReconstructedParticle, identify isolated leptons
 * based on the track cone energy, lepton identification,
 * and the track impact parameters (optional).
 */
#ifndef IsolatedPhotonFinderProcessor_h
#define IsolatedPhotonFinderProcessor_h 1

#include <string>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

class IsolatedPhotonFinderProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new IsolatedPhotonFinderProcessor ; }

		IsolatedPhotonFinderProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		/** Returns true if pfo is an isolated lepton */
		bool IsIsolatedPhoton( ReconstructedParticle* pfo ) ;

		/** Returns true if isolated, as defined by the cone energy */
		bool IsIsolatedRectangular( ReconstructedParticle* pfo ) ;
		bool IsIsolatedPolynomial( ReconstructedParticle* pfo ) ;
		bool IsIsolatedJet( ReconstructedParticle* pfo ) ;

		/** Returns true if charged */
		bool IsCharged( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes lepton ID cuts */
		bool IsPhoton( ReconstructedParticle* pfo ) ;

		/** Calculates the cone energy */
		double getConeEnergy( ReconstructedParticle* pfo ) ;

		/** [0]:Ecal energy, [1]:Hcal energy */
		void getCalEnergy( ReconstructedParticle* pfo , double* cale) ;

		/** Input collection */
		std::string _inputPFOsCollection;

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedIsoPhoCollection;

		/** Output collection of isolated leptons */
		std::string _outputIsoPhoCollection;

		LCCollection* _pfoCol;
		double _cosConeAngle;

		/** If set to true, uses PID cuts */
		bool _usePID;
		double _photonMinEcalToEcaloFraction;
		double _photonMaxEcalToEcaloFraction;

		/** If set to true, uses rectangular cuts for isolation */
		bool _useRectangularIsolation;
		double _isoMinConeEnergy;
		double _isoMaxConeEnergy;

		/** If set to true, uses polynomial cuts for isolation */
		bool _usePolynomialIsolation;
		double _isoPolynomialA;
		double _isoPolynomialB;
		double _isoPolynomialC;
		double _isoPolynomialMuonA;
		double _isoPolynomialMuonB;
		double _isoPolynomialMuonC;

		/** If set to true, uses jet-based isolation (LAL algorithm) */
		bool _useJetIsolation;
		std::string _jetCollectionName;
		std::map<ReconstructedParticle*,ReconstructedParticle*> _rpJetMap;
		double _jetIsoVetoMinXt;
		double _jetIsoVetoMaxXt;
		double _jetIsoVetoMinZ;
		double _jetIsoVetoMaxZ;
} ;

#endif

