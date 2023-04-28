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
#ifndef LGISOlatedLeptonFinderProcessor_h
#define LGISOlatedLeptonFinderProcessor_h 1

#include <string>
#include <vector>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

class LGISOlatedLeptonFinderProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new LGISOlatedLeptonFinderProcessor ; }

		LGISOlatedLeptonFinderProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		/** Returns true if pfo is an isolated lepton */
		bool IsISOlatedLepton( ReconstructedParticle* pfo ) ;

		/** Returns true if isolated, as defined by the cone energy */
		bool IsISOlatedRectangular( ReconstructedParticle* pfo ) ;
		bool IsISOlatedPolynomial( ReconstructedParticle* pfo ) ;
		bool IsISOlatedJet( ReconstructedParticle* pfo ) ;

		/** Returns true if charged */
		bool IsCharged( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes lepton ID cuts */
		bool IsLepton( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes impact parameter cuts */
		bool PassesImpactParameterCuts( ReconstructedParticle* pfo ) ; 

		/** Returns true if it passes impact parameter significance cuts */
		bool PassesImpactParameterSignificanceCuts( ReconstructedParticle* pfo ) ; 

		/** Calculates the cone energy */
		double getConeEnergy( ReconstructedParticle* pfo ) ;

		/** Get Z leptoon pair */
		vector<ReconstructedParticle*> getZLeptonPair( LCCollection* leps ) ;
		
		/** [0]:Ecal energy, [1]:Hcal energy */
		void getCalEnergy( ReconstructedParticle* pfo , double* cale) ;

		/** Input collection */
		std::string _inputPFOsCollection;

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedIsoLepCollection;
		std::string _outputPFOsRemovedIsoZLepCollection;

		/** Output collection of isolated leptons */
		std::string _outputIsoLepCollection;
		std::string _outputIsoZLepCollection;

		LCCollection* _pfoCol;
		double _cosConeAngle;

		/** If set to true, uses PID cuts */
		bool _usePID;
		double _electronMinEnergyDepositByMomentum;
		double _electronMaxEnergyDepositByMomentum;
		double _electronMinEcalToHcalFraction;
		double _electronMaxEcalToHcalFraction;
		double _muonMinEnergyDepositByMomentum;
		double _muonMaxEnergyDepositByMomentum;
		double _muonMinEcalToHcalFraction;
		double _muonMaxEcalToHcalFraction;

		/** If set to true, uses impact parameter cuts */
		bool _useImpactParameter;
		double _minD0;
		double _maxD0;
		double _minZ0;
		double _maxZ0;
		double _minR0;
		double _maxR0;
		double _mBoson;

		/** If set to true, uses impact parameter significance cuts */
		bool _useImpactParameterSignificance;
		double _minD0Sig;
		double _maxD0Sig;
		double _minZ0Sig;
		double _maxZ0Sig;
		double _minR0Sig;
		double _maxR0Sig;

		/** If set to true, uses rectangular cuts for isolation */
		bool _useRectangularIsolation;
		double _ZisoMinTrackEnergy;
		double _isoMinTrackEnergy;
		double _isoMaxTrackEnergy;
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

