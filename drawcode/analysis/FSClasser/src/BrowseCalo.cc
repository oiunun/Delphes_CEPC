#include "BrowseCalo.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <iostream>
#include <string>
#include <cmath>
#include <TVector3.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;

// protect agains rounding errors 
// will not find caps smaller than this
const float slop  = 0.25; // (mm)
const float pi    = acos(-1.0);
const float twopi = 2.0*pi;

BrowseCalo aBrowseCalo;


BrowseCalo::BrowseCalo() : Processor("BrowseCalo") {

	_description = "Performs simple digitization of sim calo hits..." ;

	std::vector<std::string> ecalCollections;

	ecalCollections.push_back(std::string("EcalBarrelCollection"));
	ecalCollections.push_back(std::string("EcalEndcapCollection"));
	ecalCollections.push_back(std::string("EcalRingCollection"));

	registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			"ECALCollections" , 
			"ECAL Collection Names" ,
			_ecalCollections ,
			ecalCollections);

	std::vector<std::string> hcalCollections;

	hcalCollections.push_back(std::string("HcalBarrelRegCollection"));
	hcalCollections.push_back(std::string("HcalEndcapRingsCollection"));
	hcalCollections.push_back(std::string("HcalEndcapsCollection"));

	registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			"HCALCollections" , 
			"HCAL Collection Names" , 
			_hcalCollections , 
			hcalCollections);

	registerOutputCollection( LCIO::CALORIMETERHIT, 
			"ECALOutputCollection" , 
			"ECAL Collection of real Hits" , 
			_outputEcalCollection , 
			std::string("ECAL")) ; 

	registerOutputCollection( LCIO::CALORIMETERHIT, 
			"HCALOutputCollection" , 
			"HCAL Collection of real Hits" , 
			_outputHcalCollection , 
			std::string("HCAL")) ; 


	registerOutputCollection( LCIO::LCRELATION, 
			"RelationOutputCollection" , 
			"CaloHit Relation Collection" , 
			_outputRelCollection , 
			std::string("RelationCaloHit")) ; 

	registerProcessorParameter("ECALThreshold" , 
			"Threshold for ECAL Hits in GeV" ,
			_thresholdEcal,
			(float)5.0e-5);

	registerProcessorParameter("HCALThreshold" , 
			"Threshold for HCAL Hits in GeV" ,
			_thresholdHcal,
			(float)2.5e-4);


	std::vector<int> ecalLayers;
	ecalLayers.push_back(20);
	ecalLayers.push_back(100);


	registerProcessorParameter("ECALLayers" , 
			"Index of ECal Layers" ,
			_ecalLayers,
			ecalLayers);



	std::vector<int> hcalLayers;
	hcalLayers.push_back(100);

	registerProcessorParameter("HCALLayers" , 
			"Index of HCal Layers" ,
			_hcalLayers,
			hcalLayers);


	std::vector<float> calibrEcal;
	calibrEcal.push_back(40.91);
	calibrEcal.push_back(81.81);


	registerProcessorParameter("CalibrECAL" , 
			"Calibration coefficients for ECAL" ,
			_calibrCoeffEcal,
			calibrEcal);


	std::vector<float> calibrHcal;
	calibrHcal.push_back(34.8);

	registerProcessorParameter("CalibrHCAL" , 
			"Calibration coefficients for HCAL" ,
			_calibrCoeffHcal,
			calibrHcal);


	registerProcessorParameter("IfDigitalEcal" ,
			"Digital Ecal" , 
			_digitalEcal , 
			0);


	registerProcessorParameter("IfDigitalHcal" ,
			"Digital Hcal" , 
			_digitalHcal , 
			0);

	registerProcessorParameter("ECALGapCorrection" , 
			"Correct for ECAL gaps" ,
			_ecalGapCorrection,
			(int)1);

	registerProcessorParameter("ECAlEndcapCorrectionFactor" , 
			"Energy correction for endcap" ,
			_ecalEndcapCorrectionFactor,
			(float)1.025);

	registerProcessorParameter("ECALGapCorrectionFactor" , 
			"Factor applied to gap correction" ,
			_ecalGapCorrectionFactor,
			(float)1.0);

	registerProcessorParameter("ECALModuleGapCorrectionFactor" , 
			"Factor applied to module gap correction" ,
			_ecalModuleGapCorrectionFactor,
			(float)0.5);

	_treeFileName="calo.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_treeName="calohits";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);


	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

}

void BrowseCalo::init() {

	_nRun = -1;
	_nEvt = 0;

	//fg: need to set default encoding in for reading old files...
	CellIDDecoder<SimCalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6") ;

   /*
	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}
	*/

	_outputHits = new TTree(_treeName.c_str(), "Calorimetry");
	_outputHits->SetAutoSave(32*1024*1024);      // autosave every 32MB

	_outputHits->Branch("Run"     ,       &_Run,            "Run/D");
	_outputHits->Branch("Evt"     ,       &_Evt,            "Evt/D");
	_outputHits->Branch("nHit"    ,      &_nHit,           "nHit/D");
	_outputHits->Branch("EnDep"   ,     &_EnDep,          "EnDep/D");
	_outputHits->Branch("EnMC"    ,      &_EnMC,           "EnMC/D");
	_outputHits->Branch("nMC"     ,       &_nMC,            "nMC/D");
	_outputHits->Branch("EnMCP"   ,      _EnMCP,      "EnMCP[99]/D");
	_outputHits->Branch("EmCali1" ,   &_EmCali1,        "EmCali1/D");
	_outputHits->Branch("EmCali2" ,   &_EmCali2,        "EmCali2/D");
	_outputHits->Branch("HdCali1" ,   &_EmCali1,        "HdCali1/D");
	_outputHits->Branch("CosTheta",   _CosTheta, "CosTheta[9999]/D");
	_outputHits->Branch("Phi"     ,        _Phi,      "Phi[9999]/D");
	_outputHits->Branch("Energy"  ,     _Energy,       "En[9999]/D");
	_outputHits->Branch("XX"      ,         _XX,       "XX[9999]/D");
	_outputHits->Branch("YY"      ,         _YY,       "YY[9999]/D");
	_outputHits->Branch("ZZ"      ,         _ZZ,       "ZZ[9999]/D");
	_outputHits->Branch("RR"      ,         _RR,       "RR[9999]/D");
	_outputHits->Branch("CellId0" ,    _CellId0,  "CellId0[9999]/D");
	_outputHits->Branch("CellId1" ,    _CellId1,  "CellId1[9999]/D");
	_outputHits->Branch("Layer"   ,      _Layer,    "Layer[9999]/D");
	_outputHits->Branch("Stave"   ,      _Stave,    "Stave[9999]/D");
	_outputHits->Branch("Module"  ,     _Module,   "Module[9999]/D");

}


void BrowseCalo::processRunHeader( LCRunHeader* run) { 

	_nRun++ ;
	_nEvt = 0;

	// Calorimeter geometry from GEAR
	const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
	const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
	const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
	const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();

	// determine geometry of ECAL
	int symmetry = pEcalBarrel.getSymmetryOrder();
	_zOfEcalEndcap = (float)pEcalEndcap.getExtent()[2];

	// Determine ECAL polygon angles
	// Store radial vectors perpendicular to stave layers in _ecalBarrelStaveDir 
	// ASSUMES Mokka Stave numbering 0 = top, then numbering increases anti-clockwise
	if(symmetry>1){
		float nFoldSymmetry = static_cast<float>(symmetry);
		float phi0 = pEcalBarrel.getPhi0();
		for(int i=0;i<symmetry;++i){
			float phi  = phi0 + i*twopi/nFoldSymmetry;
			_barrelStaveDir[i][0] = cos(phi);
			_barrelStaveDir[i][1] = sin(phi);
		}
	}  

	for(int i=0;i<ecalBarrelLayout.getNLayers();++i){
		_barrelPixelSizeT[i] = ecalBarrelLayout.getCellSize0(i);
		_barrelPixelSizeZ[i] = ecalBarrelLayout.getCellSize1(i);
	}

	for(int i=0;i<ecalEndcapLayout.getNLayers();++i){
		_endcapPixelSizeX[i] = ecalEndcapLayout.getCellSize0(i);
		_endcapPixelSizeY[i] = ecalEndcapLayout.getCellSize1(i);
	}

} 

void BrowseCalo::processEvent( LCEvent * evt ) { 

	LCCollection * col_mcp = evt->getCollection( "MCParticle" ) ;
	int _nMCP=col_mcp->getNumberOfElements();
	_EnMC=0, _nMC=0;
	for(int j = 0; j < _nMCP; j++)
	{
		MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_mcp->getElementAt(j));

		//int parentnum   = a_MCP->getParents().size();
		int daughternum = a_MCP->getDaughters().size();
		int PDGID       = a_MCP->getPDG();
		int GenStatus   = a_MCP->getGeneratorStatus();
		if( 1 == GenStatus &&  22 == abs(PDGID) ){
			_EnMCP[(int)_nMC] =  a_MCP->getEnergy();
			_EnMC        += a_MCP->getEnergy();
			_nMC++;
		}	
	}
	// create the output collections
	LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
	LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
	LCCollectionVec *relcol  = new LCCollectionVec(LCIO::LCRELATION);

	// copy the flags from the input collection
	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	ecalcol->setFlag(flag.getFlag());
	hcalcol->setFlag(flag.getFlag());

	// if making gap corrections clear the vectors holding pointers to calhits
	if(_ecalGapCorrection!=0){
		for(int is=0;is<MAX_STAVES;is++){
			for(int il=0;il<MAX_LAYERS;il++){	
				_calHitsByStaveLayer[is][il].clear();
				_calHitsByStaveLayerModule[is][il].clear();
			}
		}
	}

	_Run=evt->getEventNumber();	
	_Evt=evt->getEventNumber();	
	_EmCali1 = _calibrCoeffEcal[0];  
	_EmCali2 = _calibrCoeffEcal[1];  
	_HdCali1 = _calibrCoeffHcal[1];  
	// 
	// * Reading Collections of ECAL Simulated Hits * 
	// 
	string initString;
	int nhit = 0; 
	_EnDep   = 0;
	for (unsigned int i(0); i < _ecalCollections.size(); ++i) {
		try{
			LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
			initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
			int numElements = col->getNumberOfElements();
			CellIDDecoder<SimCalorimeterHit> idDecoder( col );
			for (int j(0); j < numElements; ++j) {
				SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
				float energy = hit->getEnergy();
				// apply threshold cut
				if (energy > _thresholdEcal) {
					CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
					int cellid  = hit->getCellID0();
					int cellid1 = hit->getCellID1();
					int layer   = idDecoder(hit)["K-1"];
					int stave   = idDecoder(hit)["S-1"];
					int module  = idDecoder(hit)["M"];
					float calibr_coeff(1.);
					// save hits by module/stave/layer if required later
					if(_ecalGapCorrection!=0){
						_calHitsByStaveLayer[stave][layer].push_back(calhit);
						_calHitsByStaveLayerModule[stave][layer].push_back(module);
					}

					// retrieve calibration constants
					for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
						int min,max;
						if (k == 0){ 
							min = 0;		      
						}else{ 
							min = _ecalLayers[k-1];
						} 
						max = _ecalLayers[k];
						if (layer >= min && layer < max) {
							calibr_coeff = _calibrCoeffEcal[k];
							break;
						}
					} 
					// apply calibration
					double cal_en = energy;
					if (_digitalEcal) {
						cal_en = calibr_coeff; 
					}
					else {
						// if in endcap apply additional factor to calibration to account for
						// the difference in response due to the orientation of B wrt absorber
						if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)energy=energy*_ecalEndcapCorrectionFactor;
						cal_en = calibr_coeff*energy;
					}
					calhit->setEnergy(cal_en); 
					// set other ECAL quanties
					calhit->setPosition(hit->getPosition());
					calhit->setType((int)0);
					calhit->setRawHit(hit);
					calhit->setCellID0(cellid);
					calhit->setCellID1(cellid1);
					ecalcol->addElement(calhit);
					// make relation between hit and sim hit
					LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.);
					relcol->addElement( rel );
            
					if( nhit > MAX_FIRED) continue;
					TVector3 v3(hit->getPosition() );
					_EnDep         += energy; //cal_en;
					_Energy[nhit]   = energy; //cal_en;
					_XX    [nhit]   = v3.X(); 
					_YY    [nhit]   = v3.Y();
					_ZZ    [nhit]   = v3.Z();
					_RR    [nhit]   = v3.Mag();
					_CellId0[nhit]  = cellid; 
					_CellId1[nhit]  = cellid1; 
					_Layer[nhit]    = layer; 
					_Stave[nhit]    = stave; 
					_Module[nhit]   = module; 
					_CosTheta[nhit] = v3.CosTheta(); 
					_Phi[nhit]      = v3.Phi(); 
					nhit++;
				}
			}
		}
		catch(DataNotAvailableException &e){ 
		}
	}

	// if requested apply gap corrections in ECAL ? 
	if(_ecalGapCorrection!=0)this->fillECALGaps();

	// add ECAL collection to event
	ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
	evt->addCollection(ecalcol,_outputEcalCollection.c_str());

	//
	// * Reading HCAL Collections of Simulated Hits * 
	//

	for (unsigned int i(0); i < _hcalCollections.size(); ++i) {
		try{
			LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
			initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
			int numElements = col->getNumberOfElements();
			CellIDDecoder<SimCalorimeterHit> idDecoder(col);
			for (int j(0); j < numElements; ++j) {
				SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
				float energy = hit->getEnergy();


				if (energy > _thresholdHcal) {

					CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
					int cellid = hit->getCellID0();
					int cellid1 = hit->getCellID1();
					float calibr_coeff(1.);
					int layer =idDecoder(hit)["K-1"]; 
					for (unsigned int k(0); k < _hcalLayers.size(); ++k) {
						int min,max;
						if (k == 0) 
							min = 0;
						else 
							min = _hcalLayers[k-1];
						max = _hcalLayers[k];
						if (layer >= min && layer < max) {
							calibr_coeff = _calibrCoeffHcal[k];
							break;
						}
					} 
					calhit->setCellID0(cellid);		  
					calhit->setCellID1(cellid1);
					double cal_en = energy;
					if (_digitalHcal) {
						cal_en = calibr_coeff; 
					}
					else {
						cal_en = calibr_coeff*energy;
					}
					calhit->setEnergy(cal_en); 
					calhit->setPosition(hit->getPosition());
					calhit->setType(int(1));
					calhit->setRawHit(hit);
					hcalcol->addElement(calhit);
					LCRelationImpl *rel = new LCRelationImpl(calhit,hit,1.0);
					relcol->addElement( rel );


					if( nhit > MAX_FIRED) continue;
					TVector3 v3(hit->getPosition() );
					_EnDep         += energy; //cal_en;
					_Energy[nhit]   = energy; //cal_en;
					_XX    [nhit]   = v3.X(); 
					_YY    [nhit]   = v3.Y();
					_ZZ    [nhit]   = v3.Z();
					_RR    [nhit]   = v3.Mag();
					_CellId0[nhit]  = cellid; 
					_CellId1[nhit]  = cellid1; 
					_Layer[nhit]    = layer; 
					_Stave[nhit]    = -1; 
					_Module[nhit]   = -1; 
					_CosTheta[nhit] = v3.CosTheta(); 
					_Phi[nhit]      = v3.Phi(); 
					nhit++;
				}

			}
		}
		catch(DataNotAvailableException &e){ 
		}
	}

	_nHit = nhit; 
	_outputHits->Fill();
	// add HCAL collection to event
	hcalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
	evt->addCollection(hcalcol,_outputHcalCollection.c_str());

	// add relation collection for ECAL/HCAL to event
	evt->addCollection(relcol,_outputRelCollection.c_str());

	_nEvt++;

}


void BrowseCalo::check( LCEvent * evt ) { }

void BrowseCalo::end(){ 


	if (_outputHits) {
		_outputHits->Write();
		//TFile *tree_file = _outputHits->GetCurrentFile(); //just in case we switched to a new file
		//tree_file->Write();
		//delete tree_file;

	}


} 


void BrowseCalo::fillECALGaps( ) { 

	// Loop over hits in the Barrel
	// For each layer calculated differences in hit positions
	// Look for gaps based on expected separation of adjacent hits
	// loop over staves and layers

	for (int is=0; is < MAX_STAVES; ++is) {
		for (int il=0; il < MAX_LAYERS; ++il) {
			if(_calHitsByStaveLayer[is][il].size()>1){
				// compare all pairs of hits just once (j>i)

				for (unsigned int i=0;i<_calHitsByStaveLayer[is][il].size()-1;++i){
					CalorimeterHitImpl* hiti = _calHitsByStaveLayer[is][il][i]; 
					int modulei = _calHitsByStaveLayerModule[is][il][i];
					float xi = hiti->getPosition()[0];
					float yi = hiti->getPosition()[1];
					float zi = hiti->getPosition()[2];

					for (unsigned int j=i+1;j<_calHitsByStaveLayer[is][il].size();++j){
						CalorimeterHitImpl* hitj = _calHitsByStaveLayer[is][il][j]; 
						int modulej = _calHitsByStaveLayerModule[is][il][j];
						float xj = hitj->getPosition()[0];
						float yj = hitj->getPosition()[1];
						float zj = hitj->getPosition()[2];
						float dz = fabs(zi-zj);
						// *** BARREL CORRECTION ***
						if( fabs(zi)<_zOfEcalEndcap && fabs(zj)<_zOfEcalEndcap){
							// account for stave directions using normals
							// calculate difference in hit postions in z and along stave
							float dx = xi-xj;
							float dy = yi-yj;
							float dt = fabs(dx*_barrelStaveDir[is][0] + dy*_barrelStaveDir[is][1]);
							// flags for evidence for gaps
							bool zgap  = false;   // in z direction
							bool tgap  = false;   // along stave 
							bool ztgap = false;  // in both z and along stave 
							bool mgap  = false;   // gaps between ECAL modules

							// criteria gaps in the z and t direction
							float zminm = 1.0*_barrelPixelSizeZ[il]-slop;
							float zmin  = 1.0*_barrelPixelSizeZ[il]+slop;
							float zmax  = 2.0*_barrelPixelSizeZ[il]-slop;
							float tminm = 1.0*_barrelPixelSizeT[il]-slop;
							float tmin  = 1.0*_barrelPixelSizeT[il]+slop;
							float tmax  = 2.0*_barrelPixelSizeT[il]-slop;

							// criteria for gaps
							// WOULD BE BETTER TO USE GEAR TO CHECK GAPS ARE OF EXPECTED SIZE
							if( dz > zmin  && dz < zmax && dt < tminm )zgap = true;
							if( dz < zminm && dt > tmin && dt < tmax )tgap = true;
							if( dz > zmin  && dz < zmax && dt > tmin && dt < tmax )ztgap=true;

							if(modulei!=modulej){
								if( dz > zmin && dz < 3.0*_barrelPixelSizeZ[il]-slop && dt < tmin)mgap = true;
							}

							// found a gap now apply a correction based on area of gap/area of pixel
							if(zgap||tgap||ztgap||mgap){
								float ecor = 1.;
								float f = _ecalGapCorrectionFactor; // fudge
								if(mgap)f = _ecalModuleGapCorrectionFactor;
								if(zgap||mgap)ecor = 1.+f*(dz - _barrelPixelSizeZ[il])/2./_barrelPixelSizeZ[il];
								if(tgap)ecor = 1.+f*(dt - _barrelPixelSizeT[il])/2./_barrelPixelSizeT[il];
								if(ztgap)ecor= 1.+f*(dt - _barrelPixelSizeT[il])*(dz - _barrelPixelSizeZ[il])/4./_barrelPixelSizeT[il]/_barrelPixelSizeZ[il];     
								float ei = hiti->getEnergy()*ecor;
								float ej = hitj->getEnergy()*ecor;
								hiti->setEnergy(ei);
								hitj->setEnergy(ej);
							}

							// *** ENDCAP CORRECTION ***
						}else if(fabs(zi)>_zOfEcalEndcap && fabs(zj)>_zOfEcalEndcap&&dz<100){
							float dx = fabs(xi-xj);
							float dy = fabs(yi-yj);
							bool xgap = false;
							bool ygap = false;
							bool xygap = false;
							// criteria gaps in the z and t direction
							float xmin  = 1.0*_endcapPixelSizeX[il]+slop;
							float xminm = 1.0*_endcapPixelSizeX[il]-slop;
							float xmax  = 2.0*_endcapPixelSizeX[il]-slop;
							float ymin  = 1.0*_endcapPixelSizeY[il]+slop;
							float yminm = 1.0*_endcapPixelSizeY[il]-slop;
							float ymax  = 2.0*_endcapPixelSizeY[il]-slop;
							// look for gaps
							if(dx > xmin && dx < xmax && dy < yminm )xgap = true;
							if(dx < xminm && dy > ymin && dy < ymax )ygap = true;
							if(dx > xmin && dx < xmax && dy > ymin && dy < ymax )xygap=true;

							if(xgap||ygap||xygap){
								// found a gap make correction
								float ecor = 1.;
								float f = _ecalGapCorrectionFactor; // fudge
								if(xgap)ecor = 1.+f*(dx - _endcapPixelSizeX[il])/2./_endcapPixelSizeX[il];
								if(ygap)ecor = 1.+f*(dy - _endcapPixelSizeY[il])/2./_endcapPixelSizeY[il];
								if(xygap)ecor= 1.+f*(dx - _endcapPixelSizeX[il])*(dy - _endcapPixelSizeY[il])/4./_endcapPixelSizeX[il]/_endcapPixelSizeY[il];     
								hiti->setEnergy( hiti->getEnergy()*ecor );
								hitj->setEnergy( hitj->getEnergy()*ecor );
							}
						}
					}
				}
			}
		}
	}

	return;

}

