#ifndef DIGITIZER_NEWBROWSECALO_H
#define DIGITIZER_NEWBROWSECALO_H 1

#include "marlin/Processor.h"
#include <IMPL/CalorimeterHitImpl.h>
#include <marlinutil/CalorimeterHitType.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>

using namespace lcio ;
using namespace marlin ;

const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;
const int MAX_FIRED  = 99999;

class NewBrowseCalo : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new NewBrowseCalo ; }
  
  
  NewBrowseCalo() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
   
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

  virtual void fillECALGaps() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _ecalCollections;
  std::vector<std::string> _hcalCollections;


  std::string _outputEcalCollection0;
  std::string _outputEcalCollection1;
  std::string _outputEcalCollection2;
  std::string _outputHcalCollection0;
  std::string _outputHcalCollection1;
  std::string _outputHcalCollection2;
  std::vector<std::string> _outputEcalCollections;
  std::vector<std::string> _outputHcalCollections;
  std::string _outputRelCollection;

  float _thresholdEcal;
  std::vector<float> _thresholdHcal;

  int _digitalEcal;
  int _digitalHcal;

  std::vector<float> _calibrCoeffEcal;
  std::vector<float> _calibrCoeffHcal;

  std::vector<int> _ecalLayers;
  std::vector<int> _hcalLayers;

  int _ecalGapCorrection;
  float _ecalGapCorrectionFactor;
  float _ecalModuleGapCorrectionFactor;
  float _ecalEndcapCorrectionFactor;
  float _hcalEndcapCorrectionFactor;

  std::vector<CalorimeterHitImpl*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap;
  float _barrelPixelSizeT[MAX_LAYERS];
  float _barrelPixelSizeZ[MAX_LAYERS];
  float _endcapPixelSizeX[MAX_LAYERS];
  float _endcapPixelSizeY[MAX_LAYERS];
  float _barrelStaveDir[MAX_STAVES][2];



  std::string _treeFileName;
  std::string _treeName;

  int    _overwrite;
  TTree *_outputHits; 

  std::ostream *_output;
  std::string   _fileName;
  std::string   _histFileName;

  double _Run, _Evt, _nHit, _EnDep, _EnMC, _nMC;
  int    _Type[MAX_FIRED];
  double _EmCali1, _EmCali2, _HdCali1;
  double _EnMCP[100];
  double _XX[MAX_FIRED],        _YY[MAX_FIRED], _ZZ[MAX_FIRED];
  double _CosTheta[MAX_FIRED], _Phi[MAX_FIRED], _RR[MAX_FIRED];
  double _Energy[MAX_FIRED];
  double _CellId0[MAX_FIRED];
  double _CellId1[MAX_FIRED];
  double _Layer  [MAX_FIRED];
  double _Stave  [MAX_FIRED];
  double _Module [MAX_FIRED];



} ;

#endif



