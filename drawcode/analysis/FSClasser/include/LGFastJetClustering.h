#ifndef LGFastJetClustering_h
#define LGFastJetClustering_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include "IMPL/LCFlagImpl.h"

#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include <string>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include "TLorentzVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Utility.h"
#include "cnpy.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "EnergyCorrelator.hh" // In external code, this should be fastjet/contrib/EnergyCorrelator.hh

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace Utility;
using namespace EVENT;
//using namespace fastjet;
using namespace fastjet::contrib;

class LGFastJetClustering : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new LGFastJetClustering ; }
    
  LGFastJetClustering() ;
  
  virtual void init() ;

  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void processJets ( LCEvent * evt , LCCollectionVec* jets ) ; 
    
  virtual void check( LCEvent * evt ) ; 

  virtual void end() ;

  fastjet::JetAlgorithm  GetAlgorithm()  const {return fAlgorithm;}

  void SetAlgorithm(fastjet::JetAlgorithm f)  {fAlgorithm = f;}

  void thrust(const vector<fastjet::PseudoJet> sortedInputs ) ; 
  void fox_wolfram(const vector<fastjet::PseudoJet> sortedInputs ) ; 
  MCParticle * FindParton(MCParticle *mcp); 
  double * angle (const TLorentzVector& bbst)
  {
	  static double ang[2]={0,0};
	  ang[0]=0; ang[1]=0;
	  double rxy=sqrt(pow(bbst.X(),2.)+pow(bbst.Y(),2.));
	  if(bbst.Rho()<=1e-6)
	  {
	  }
	  else if(rxy<=1e-6)
	  {
		  if(bbst.Z()<0.0) ang[0]=3.1415927;
	  }
	  else
	  {
		  double theta=acos(bbst.Z()/bbst.Rho());
		  double phi = (bbst.X() == 0.0 && bbst.Y() == 0.0 && bbst.Z() == 0.0 ? 0.0:TMath::ATan2(bbst.Y(),bbst.X()));
		  ang[0]=theta;
		  ang[1]=phi;
	  }
	  return (double*)ang;
  }

  TLorentzVector erout4(const TLorentzVector& p1, double phi, double theta )
  {
	  double cp = cos(phi);
	  double sp = sin(phi);
	  double ct = cos(theta);
	  double st = sin(theta);
	  double t  = p1.E(), x=p1.X(), y=p1.Y(), z=p1.Z();
	  double xp=x*cp*ct+y*sp*ct-z*st, yp= -x*sp+y*cp, zp=x*cp*st+y*sp*st+z*ct;
	  TLorentzVector Rp(zp,xp,yp,t);
	  return Rp;
  }

  void sort_indices(vector<int> & indices, 
		  const vector<double> & values) {
	  fastjet::IndexedSortHelper index_sort_helper(&values);
	  sort(indices.begin(), indices.end(), index_sort_helper);
  }

  void analyze(const vector<fastjet::PseudoJet> & Jets ) ;

  template<class T> vector<T>  objects_sorted_by_values(
		  const vector<T> & objects, 
		  const vector<double> & values) {

	  assert(objects.size() == values.size());

	  // get a vector of indices
	  vector<int> indices(values.size());
	  for (size_t i = 0; i < indices.size(); i++) {indices[i] = i;}

	  // sort the indices
	  sort_indices(indices, values);

	  // copy the objects 
	  vector<T> objects_sorted(objects.size());

	  // place the objects in the correct order
	  for (size_t i = 0; i < indices.size(); i++) {
		  objects_sorted[i] = objects[indices[i]];
	  }

	  return objects_sorted;
  }

  //----------------------------------------------------------------------
  /// return a vector of jets sorted into decreasing energy
  vector<ReconstructedParticle*> sorted_by_E(const vector<ReconstructedParticle*> & jets) {
	  vector<double> energies(jets.size());
	  for (size_t i = 0; i < jets.size(); i++) {energies[i] = -jets[i]->getEnergy();}
	  return objects_sorted_by_values(jets, energies);
  }



 protected:

  // root file and tree objects
  TFile * _rootfile;
  TTree * _Etree;
  TTree * _Ctree;
  TTree * _Ftree;

  std::vector<double> dat;
  std::vector<double> tag ;
  long unsigned int Ne, Nv, Pm; 

  LCFlagImpl Cluflag;

  std::string _inputCollection;
  std::string _inputMCTruthMap;
  std::string _outputCollection;
  std::string _outputIsoLepCol;
  std::string _outputRemainCol;
  
  LCCollectionVec* _jetsCol, *_leptCol, *_leftCol;

  fastjet::JetAlgorithm fAlgorithm; 

  std::string sAlgorithm; 

  // px, py, pz, E, nptc
  float _jetVector[5][50];

  float _eCMS;

  int _nRun, _nEvt, _nJets, _nPFOmin, _nJetsHE, _BCL, _NPJ, _label, _PDGID[200], _Charge[200];
  int _lab_ll, _lab_cc, _lab_bb, _lab_ee, _lab_uu, _lab_tt, _lab_gg, _lab_aa, _lab_zz, _lab_ww, _lab_az; 
  double  _THETA[200],  _PHI[200], _CosT[200], _SinT[200], _delR[200], _PhiSinT[200];
  double  _Energy[200], _Mom[200], _logEn[200], _logMo[200],_logEnF[200], _logMoF[200];
  double  _D0[200], _Z0[200];
  double  _llRec, _llMass, _hMass, _hRecMass, _TrueHmass, _Hmass,  _TrueZmass, _Zmass, _TrueWmass, _Wmass; 
  double  _Sphericity ;
  double  _Aplanarity ;
  double  _C          ;
  double  _D          ;
  double  _Thrust, _Theta, _Phi, _ThrEDM;
  double  _Major,  _Minor, _Thetap, _Phip, _ThrEDMp;  
  double  _FM [50];  
  double  _FMI[50];  
  double  _FMS[50];  
  double  _Ci[28];  
  double  _JetCharge, _JetPDGID;

  int _print, _nJetMax, _fillTree;

  vector<string> _RemoveList;
  vector<int>    _PdgidList;
  int            _Charged, _InclusiveExclusive, _save_npz;
  double RemoveResonance, _EnergyCut; 

  double _RPar, _rp, _eJet, _PPar, _pp, _PtCut;
} ;

#endif
