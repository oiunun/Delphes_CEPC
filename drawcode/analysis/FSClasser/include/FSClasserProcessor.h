#ifndef FSClasserProcessor_h
#define FSClasserProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "CLHEP/Units/PhysicalConstants.h"
//#include "CLHEP/Geometry/Point3D.h"
//#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Vector/LorentzVector.h"
//#include "CLHEP/Vector/TwoVector.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include <TVectorD.h> 
#include <TMatrixDSym.h> 
#include <TMatrixDSymEigen.h> 
#include <TTree.h>
#include <MCTruthHelper.h>
#include <NTupleHelper.h>
#include <FSHelper.h>

#include "BaseFitObject.h"
#include "ParticleFitObject.h"
#include "JetFitObject.h"
#include "ISRPhotonFitObject.h"
#include "PConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "FourJetPairing.h"
#include "MassConstraint.h"
#include "cepcplotstyle.h"
#include "Utilities.h"



using namespace lcio ;
using namespace marlin ;

class TTree;
class MCTruthHelper;
class NTupleHelper;

/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: FSClasserProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class FSClasserProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FSClasserProcessor ; }
  
  
  FSClasserProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
 
  void    BuildFullSimParticleList(
		  LCCollection    *colMC  , LCCollection    *col2Jet, LCCollection    *col4Jet, 
		  LCCollection   *colMCTL , LCCollection    *colPFO ,  
		  LCCollection  *colPFOPandora=NULL, LCCollection  *colIsoLep=NULL,
      LCCollection  *colTaus=NULL
		  ); 
  bool    IsPhoton( ReconstructedParticle *par, LCCollection *pfoCol ); 
  int     IsLepton( ReconstructedParticle *pfo ); 
  void    getConeNeutraChargeInfo(ReconstructedParticle *part, LCCollection *collPfo, 
		  double cosCone, int &nConeCharged, int &nConeNeutral,
		  double &coneEC, double &coneEN);
  void    getCalEnergy(ReconstructedParticle *part, double *cale);
  double  getTransverseMomentum(ReconstructedParticle *part);
  void    CleanEvt() ; 
  void    MakePlots() ; 
  bool    checkCombination  (const vector<FSParticle*>& combo, bool complete, bool inclusive);
  double  getMassColApp(TLorentzVector p4plus, TLorentzVector p4minus, double missingx, double missingy  );
  double  sphericity(); 
  double  Lsphericity(); 
  double  thrust(); 
  double  FDevent();
  int     NHScale( std::vector<TVector3>, double, double);	

 protected:

  static const unsigned int  MAXFS  = 199;
  /** Input collection name.
   */
  MCTruthHelper* m_mcTruthHelper;
  NTupleHelper*  m_ntGEN;
  NTupleHelper*  m_ntANA;
  NTupleHelper*  m_ntVal;

  std::string   _treeFileName;
  std::string   _treeName;

  std::string   _colMCP  ;
  std::string   _colPFOs ;
  std::string   _colNewPFOs ;
  std::string   _colIsoLeps ;
  std::string   _colIsoZLeps;
  std::string   _colTaus ;
  std::string   _col2Jets;
  std::string   _col4Jets;
  std::string   _colLeptons ;
  std::string   _colMCTL ;
  TH1D *hValid[30];
  TH1D *hPando[30];
  TH1D *hRes  [ 2];

  int _nRun ;
  int _nEvt ;
  int _nb, _nc, _nl;
  int _nhfs, _nhpp, _nhww, _nhzz, _nhee, _nhmm, _nhtt, _nhuu, _nhdd, _nhss, _nhbb, _nhcc, _nhgg;

  TLorentzVector pVis, pVis_Pandora, totalP4;
  TLorentzVector m_ecms;
  int ntrks_Pandora, nclus_Pandora, nPFOs_Pandora;
  int ntrks        , nclus        , nPFOs        , numberJets , numberTaus;

  double _y12, _y23, _y34, _y45, _y56, _y67, _Emax, _Pmax, _FD;
  double b, ISRPzMaxB; 

  //int          _overwrite;
  //TFile        *tree_file;
  //TTree        *_outputTree;

  vector<FSInfo*>   m_FSInfoVector;
  int     m_debug ;
  int     m_savemc;
  int     m_makeplots;
  int     m_matchmc;
  int     m_kmfit ;
  int     m_TagFlavor;
  int     m_showmc;
  int     m_luxury;
  int     m_savehis;
  int     m_full;
  int     m_useCalo;
  int     m_EventShape, m_LinearSphericity;
  double  m_ECM;
  double  m_wate;
  double  m_kappa;
  double  m_EnergyCut, m_LepEnergyThreshold ;
  double  m_Sphericity, m_Aplanarity, m_C, m_D, m_Thrust, m_Theta, m_Phi, m_ThrEDM;
  double  m_Major, m_Minor, m_Thetap, m_Phip, m_ThrEDMp;  
  bool    m_findEp;
  bool    m_findEm;
  bool    m_findMup;
  bool    m_findMum;
  bool    m_findTaup;
  bool    m_findTaum;
  bool    m_findJets;
  bool    m_findPhoton;
  
  vector<FSParticle*>    mupList      ;
  vector<FSParticle*>    mumList      ;
  vector<FSParticle*>     epList      ;
  vector<FSParticle*>     emList      ;
  vector<FSParticle*>   taupList      ;
  vector<FSParticle*>   taumList      ;
  vector<FSParticle*>    jetList      ;
  vector<FSParticle*>   jet4List      ;
  vector<FSParticle*> photonList      ;
  vector<FSParticle*> ParticleTrash   ;
  map< string, vector<FSParticle*> > DictPList;
  vector<ReconstructedParticle*>  raw_PFOs; 
  vector<Track*>  raw_trackvec;
  TVector3         m_HiggsBoostVector; 

  
  int               m_CutPass [MAXFS];
  double            m_Checking[MAXFS];
  string            m_FSNames [MAXFS];
  vector<string>    m_FSCuts  [MAXFS];
  string            m_FSMaCs  [MAXFS];
  string            m_FSMMFits[MAXFS];
  double            m_cutpass [MAXFS][20];
  string            m_boostFS;

} ;

#endif



