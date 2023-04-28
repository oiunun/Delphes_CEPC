/**
 */
#ifndef ISRFinderProcessor_h
#define ISRFinderProcessor_h 1

#include <string>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;

class ISRFinderProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new ISRFinderProcessor ; }

		ISRFinderProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;

	protected:

		/** Returns true if pfo is an initial state radiation */
		bool IsISR( ReconstructedParticle* pfo ) ;

		/** Returns true if neutral */
		bool IsCharged( ReconstructedParticle* pfo ) ;

		/** Returns true if it passes energy and cos theta cuts */
		bool IsEnergyandCos( ReconstructedParticle* pfo ) ;

		/** Retruns true if it deposit most of its energy in ECAL*/
		bool PhotonSelection( ReconstructedParticle* pfo ) ;

		/** Input collection */
		std::string _inputPFOsCollection;

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedISRCollection;

		/** Output collection of isolated leptons */
		std::string _outputISRCollection;

		LCCollection* _pfoCol;
} ;

#endif

