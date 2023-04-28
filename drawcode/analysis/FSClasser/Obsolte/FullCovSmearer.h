#ifndef FullCovSmearer_h
#define FullCovSmearer_h 1

#include "marlin/MarlinConfig.h"

#ifdef MARLIN_CLHEP  // only if CLHEP is available !

#include "CLHEP/Vector/LorentzVector.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/MCParticleImpl.h"

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace CLHEP ;
using namespace IMPL ;
using namespace EVENT ;


namespace marlin{


  /** Interface for smearing of four vectors - based on CLHEP::HepLorentzVector
   *
   *  @author F. Gaede, DESY
   *  @version $Id: FullCovSmearer.h,v 1.3 2006-03-30 16:12:16 gaede Exp $ 
   */ 
  
  class FullCovSmearer {
    
  public:
    
    /** Virtual d'tor.*/
    virtual ~FullCovSmearer() {} 
    
    /** Smears the given four vector 
     */ 
    virtual HepLorentzVector smearedFourVector( const HepLorentzVector&  v, const int pdgCode )  =0 ;
    virtual void    smearedTrack     ( const MCParticle * mcp, const int pdgCode, TrackImpl   & trk )  = 0;
    virtual void    smearedCluster   ( const MCParticle * mcp, const int pdgCode, ClusterImpl & clu )  = 0;
    
  } ;
  
  
} // end namespace 

#endif // MARLIN_CLHEP
#endif // FullCovSmearer_h

