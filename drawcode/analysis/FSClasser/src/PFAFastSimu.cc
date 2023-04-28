/*

	A PFA based Fast Simulation tool.

Input: MCParticle - Final State MCParticles
Output: Reconstructed Particle (Arbor Objects)

Logic: For all Final State MCParticles, considering the 
Reconstruction efficiency
and 
Energy/Momentum Smearing
as well as
Confusion

according to the 

Particle Type: Neutrino, Electron-Positron, Charged Hadrond & Muon, Neutral Hadrons
and
Energy/Polar angle Distribution


For Charged Particle, 

Chance of Splitting is taken into account.

For Neutral Particle, calculate the minimal distance to a NEARBY Charged Hadron-Muon

Change of Merge is taken into account.


Author: Manqi

Nov. 28, 2013

*/



#include <PFAFastSimu.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ClusterImpl.h"
#include "IMPL/TrackImpl.h"
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h> 
#include <fstream>		
#include <cmath>
#include <vector>
#include <TF1.h>
#include <TMath.h>

using namespace std;
//PDG: gluon = 21, gamma = 22, Z = 23, W = 24, H = 25

const double mZ = 91.2;
const double mW = 80.4; 
const double mH = 125; 
const double JetReso = 0.04;	// 4% of jet energy resolution
const double PhotoReso = 0.20; // 0.16; 
const double PhotoResoConst = 0.01; 
const double NHReso = 0.60;
const double Zwidth = 2.495;
const double Wwidth = 2.085;
const double Hwidth = 0.004;	//SM Higgs Width
const double sqrtS = 250.0;
double Response[5] = {0, 0, 0, 0, 0};


double ealpha[10];
double en[10];
double emean[10];
double esigma[10];

double mumean1[10][19];
double mumean2[10][19];
double musigma1[10][19];
double musigma2[10][19];
double mufrac[10][19];
double mueff[10][19];


double pimean1[10][19];
double pimean2[10][19];
double pisigma1[10][19];
double pisigma2[10][19];
double pifrac[10][19];
double pieff[10][19];


double CrystalBall(double* x, double* par)
{
	double alpha = par[0];
	double n = par[1];
	double mean = par[2];
	double sigma = par[3];
	double A = TMath::Power(n/TMath::Abs(alpha), n)*TMath::Exp(-alpha*alpha/2.);
	double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);
	double fval = 0.0;
	double tmp = (x[0]-mean)/sigma;
	if(tmp>-alpha) {
		fval = TMath::Exp(-tmp*tmp/2.);
	}
	else {
		fval = A*TMath::Power(B-tmp, -n);
	}
	return fval;
}

PFAFastSimu a_PFAFastSimu_instance;

PFAFastSimu::PFAFastSimu()
	: Processor("PFAFastSimu"),
	_output(0)
{
	_description = "Print MC Truth" ;


	registerInputCollection( LCIO::MCPARTICLE,
			"InputCollectionName" , 
			"Name of the MCParticle input collection"  ,
			_inputCollectionName ,
			std::string("MCParticle") ) ;


	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			"RecoParticleCollectionName" , 
			"Name of the ReconstructedParticles output collection"  ,
			_recoParticleCollectionName ,
			std::string("ReconstructedParticles") ) ;

	registerOutputCollection( LCIO::LCRELATION,
			"MCTruthMappingCollectionName" , 
			"Name of the MCTruthMapping output collection"  ,
			_mcTruthCollectionName ,
			std::string("MCTruthMapping") ) ;


	_treeFileName="MCTruth.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_treeName="MCPart";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);

	_EnableSmearing = 1;
	registerProcessorParameter( "EnableSmearing" ,
			"Enable Smearing to Final State Particle Momenta." ,
			_EnableSmearing ,
			_EnableSmearing);

	_EnableEfficiency = 1;
	registerProcessorParameter( "EnableEfficiency" ,
			"Enable Geometry/Energy Dependent Efficiency." ,
			_EnableEfficiency ,
			_EnableEfficiency);

	_EnableChargeSplitting = 1; 
	registerProcessorParameter( "EnableChargeSplitting" ,
			"Enable PFA double counting" ,
			_EnableChargeSplitting ,
			_EnableChargeSplitting);

	_EnableNeutralMerge = 1; 
	registerProcessorParameter( "EnableNeutralMerge" ,
			"Enable absobtion of Neutral Cluster to Charged." ,
			_EnableNeutralMerge ,
			_EnableNeutralMerge);

	_EnableThetaScan = 1; 
	registerProcessorParameter( "EnableThetaScan" ,
			"Enable scan the polar angle" ,
			_EnableThetaScan ,
			_EnableThetaScan);

	_overwrite=1;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

	_flavor=0;
	registerProcessorParameter( "Flavor" , 
			"quark flavor to be selected" ,
			_flavor ,
			_flavor);

	_rejectNeutrino=1;
	registerProcessorParameter( "RejectNeutrino" , 
			"reject the undetectabele final state particles" ,
			_rejectNeutrino ,
			_rejectNeutrino);
	
	_momentumCut=0.3;
	registerProcessorParameter( "MomentumCut" , 
			"reject the very low energy paticles" ,
			_momentumCut ,
			_momentumCut);
}

double InvMass(double* P)
{
	double mass = 0; 

	mass = sqrt(fabs(P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - P[3]*P[3] ));

	return mass; 
}

double* DetectorReso(MCParticle * a_MCP)	//Relative scale
{

	double Reso = 0;
	double Efficiency = 1.0;
	double SplitChance = 0;

	TVector3 currP = a_MCP->getMomentum();
	int currPID = a_MCP->getPDG();
	double PAmplitude = currP.Mag();
	double cosTheta = fabs(a_MCP->getMomentum()[2]/PAmplitude);
	double CBAlpha=0;
	double CBN=0;
	double kResoLow=0;
	double kResoHigh=0;
	double kAlphaLow=0;
	double kAlphaHigh=0;
	double kNLow=0;
	double kNHigh=0;
	double ObjkEff = 1.;
	double ObjkReso = 0;
	double ObjkAlpha = 0;
	double ObjkN = 0;
	double EnDisToLow = 0;
	double EnDisToHigh = 0;
	double CoDisToLow = 0;
	double CoDisToHigh = 0;
	double Scalefactor = 0;  //Correct EndCap divergence of Charged Particle...

	double RefEnergy[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	double RefCostheta[19] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00};
	//        double ekCoeff[10] = {0.0145464, 0.0386866, 0.0663844, 0.0984627, 0.139471, 0.179123, 0.219332, 0.276281, 0.331001, 0.382948};  // for electron Radius = 1365 mm
	//        double mkCoeff[10] = {0.0153967, 4.08206e-02, 6.91961e-02, 1.03339e-01, 1.41198e-01, 1.82501e-01, 2.26472e-01, 2.82907e-01, 3.36424e-01, 3.99089e-01}; // for muon Radius = 1365 mm

	if( abs(currPID) == 12 || abs(currPID) == 14 || abs(currPID) == 16 )
	{
		Efficiency = 0;                       //Only a flat
	}
	else if(abs(currPID) == 11 )
	{
		if(PAmplitude < 10)
		{
			ObjkReso = 0.015099;
			ObjkAlpha = 1.3659;
			ObjkN = 0.99108;
		}
		else if(PAmplitude > 100)
		{
			ObjkReso = 0.39457;
			ObjkAlpha = 1.5207;
			ObjkN = 0.93348;
		}
		else
		{
			for(int i = 0; i < 9; i++)
			{
				if(PAmplitude < RefEnergy[i+1] && PAmplitude > RefEnergy[i])
				{
					kResoLow = esigma[i];
					kResoHigh = esigma[i+1];
					kAlphaLow = ealpha[i];
					kAlphaHigh = ealpha[i+1];
					kNLow = en[i];
					kNHigh = en[i+1];
					EnDisToLow = PAmplitude - RefEnergy[i];
					EnDisToHigh = RefEnergy[i+1] - PAmplitude;
					ObjkReso = (kResoLow*EnDisToHigh + kResoHigh*EnDisToLow)/(RefEnergy[i+1] - RefEnergy[i]);
					ObjkAlpha = (kAlphaLow*EnDisToHigh + kAlphaHigh*EnDisToLow)/(RefEnergy[i+1] - RefEnergy[i]);
					ObjkN = (kNLow*EnDisToHigh + kNHigh*EnDisToLow)/(RefEnergy[i+1] - RefEnergy[i]);
				}
			}
		}

		if(cosTheta > 0.86)     //Scale as effective R^2
		{
			Scalefactor = (1.0/(cosTheta*cosTheta) - 1)*2.96; //2.96 = (Half_Z/Radius)**2
			ObjkReso = ObjkReso*1.0/Scalefactor;
		}

		Reso = ObjkReso;         // Relative
		CBAlpha = ObjkAlpha;
		CBN = ObjkN;
	}
	else if(abs(currPID) == 13)
	{
		if(PAmplitude < 10)
		{
			ObjkReso = (gRandom->Gaus(PAmplitude, sqrt(mufrac[0][9]*musigma1[0][9]*musigma1[0][9]+(1-mufrac[0][9])*musigma2[0][9]*musigma2[0][9])))-PAmplitude;
			ObjkEff = mueff[0][9];
		}
		else if(PAmplitude > 100)
		{
			ObjkReso = (gRandom->Gaus(PAmplitude, sqrt(mufrac[0][9]*musigma1[0][9]*musigma1[0][9]+(1-mufrac[0][9])*musigma2[0][9]*musigma2[0][9])))-PAmplitude;
			ObjkEff = mueff[9][9];
		}
		else
		{
			for(int i = 0; i < 9; i++)
			{
				if(PAmplitude < RefEnergy[i+1] && PAmplitude > RefEnergy[i])
				{
					EnDisToLow = PAmplitude - RefEnergy[i];
					EnDisToHigh = RefEnergy[i+1] - PAmplitude;
					int id = i;
					if(EnDisToLow>EnDisToHigh) id=i+1;
					if(cosTheta<0.1) 
					{
						if(gRandom->Rndm(1)<mufrac[id][0])
							ObjkReso = (gRandom->Gaus(PAmplitude, musigma1[id][0]))-PAmplitude;
						else ObjkReso = (gRandom->Gaus(PAmplitude, musigma2[id][0]))-PAmplitude;
						ObjkEff = mueff[id][0];
					}
					else
					{
						for(int j = 0; j < 19; j++)
						{
							if(cosTheta < RefCostheta[j+1] && cosTheta > RefCostheta[j])
							{
								CoDisToLow = cosTheta - RefCostheta[i];
								CoDisToHigh = RefCostheta[i+1] - cosTheta;
								int ID = j;
								if(CoDisToLow>CoDisToHigh) ID=j+1;
								if(gRandom->Rndm(1)<mufrac[id][ID])
									ObjkReso = (gRandom->Gaus(PAmplitude, musigma1[id][ID]))-PAmplitude;
								else ObjkReso = (gRandom->Gaus(PAmplitude, musigma2[id][ID]))-PAmplitude;
								ObjkEff = mueff[id][ID];
							}
						}
					}
				}
			}
		}

		Reso = ObjkReso;         // Relative
		Efficiency = ObjkEff;
	}
	else if(abs(currPID) == 211)
	{
		if(PAmplitude < 10)
		{
			ObjkReso = (gRandom->Gaus(PAmplitude, sqrt(pifrac[0][9]*pisigma1[0][9]*pisigma1[0][9]+(1-pifrac[0][9])*pisigma2[0][9]*pisigma2[0][9])))-PAmplitude;
			ObjkEff = pieff[0][9];
		}
		else if(PAmplitude > 100)
		{
			ObjkReso = (gRandom->Gaus(PAmplitude, sqrt(pifrac[0][9]*pisigma1[0][9]*pisigma1[0][9]+(1-pifrac[0][9])*pisigma2[0][9]*pisigma2[0][9])))-PAmplitude;
			ObjkEff = pieff[9][9];
		}
		else
		{
			for(int i = 0; i < 9; i++)
			{
				if(PAmplitude < RefEnergy[i+1] && PAmplitude > RefEnergy[i])
				{
					EnDisToLow = PAmplitude - RefEnergy[i];
					EnDisToHigh = RefEnergy[i+1] - PAmplitude;
					int id = i;
					if(EnDisToLow>EnDisToHigh) id = i+1;
					if(cosTheta<0.1) 
					{
						if(gRandom->Rndm(1)<pifrac[id][0])
							ObjkReso = (gRandom->Gaus(PAmplitude, pisigma1[id][0]))-PAmplitude;
						else ObjkReso = (gRandom->Gaus(PAmplitude, pisigma2[id][0]))-PAmplitude;
						ObjkEff = pieff[id][0];
					}
					else
					{
						for(int j = 0; j < 19; j++)
						{
							if(cosTheta < RefCostheta[j+1] && cosTheta > RefCostheta[j])
							{
								CoDisToLow = cosTheta - RefCostheta[i];
								CoDisToHigh = RefCostheta[i+1] - cosTheta;
								int ID = j;
								if(CoDisToLow>CoDisToHigh) ID = j+1;
								if(gRandom->Rndm(1)<pifrac[id][ID])
									ObjkReso = (gRandom->Gaus(PAmplitude, pisigma1[id][ID]))-PAmplitude;
								else ObjkReso = (gRandom->Gaus(PAmplitude, pisigma2[id][ID]))-PAmplitude;
								ObjkEff = pieff[id][ID];
							}
						}
					}
				}
			}
		}

		Reso = ObjkReso;         // Relative
		Efficiency = ObjkEff;
	}
	else if(currPID == 22)
	{
		Reso = (gRandom->Gaus(PAmplitude,sqrt(PhotoReso*PhotoReso/PAmplitude + PhotoResoConst*PhotoResoConst)))-PAmplitude;
	}
	else	//Neutral Hadron - polar dependence?
	{
		Reso = sqrt(NHReso*NHReso/PAmplitude);                    //Protection
	}

	Response[0] = Reso;
	Response[1] = Efficiency; 
	Response[2] = SplitChance;
	Response[3] = CBAlpha;
	Response[4] = CBN; 

	return Response;
}

void PFAFastSimu::init() {

	printParameters();
	
	_nRun = 0 ;
	_nEvt = 0 ;
	
	char* Path, file[100];
	Path = getenv ("FSClasser_HOME");
   
   sprintf(file,"%s/%s", Path,"parameter/electron.para");
	ifstream electron(file);
	for(int i=0; i<10; i++)
	{
		electron>>ealpha[i]>>en[i]>>emean[i]>>esigma[i];
	}
	electron.close();

   sprintf(file,"%s/%s", Path,"parameter/result.txt");
	ifstream muon(file);
	for(int i=0; i<10; i++)
	{
		for(int j=0; j<19; j++)
		{
			muon>>mufrac[i][j]>>mumean1[i][j]>>musigma1[i][j]>>mumean2[i][j]>>musigma2[i][j]>>mueff[i][j];
		}
	}
	muon.close();


   sprintf(file,"%s/%s", Path,"parameter/pion.para");
	ifstream pion(file);
	for(int i=0; i<10; i++)
	{
		for(int j=0; j<19; j++)
			pion>>pifrac[i][j]>>pimean1[i][j]>>pisigma1[i][j]>>pimean2[i][j]>>pisigma2[i][j]>>pieff[i][j];
	}
	pion.close();

}

void PFAFastSimu::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{		
		try 	
		{    
			FreeDelAll(_tracktrash);
			FreeDelAll(_clustertrash);
			FreeDelAll(_particletrash);

			_eventNr=evtP->getEventNumber();

			//std::cout<<"Num "<<_nEvt<<std::endl;

			LCCollectionVec *fastrecoparticle = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

			LCRelationNavigator relNav( LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;

			Cluflag.setBit(LCIO::CHBIT_LONG);

			fastrecoparticle->setFlag(Cluflag.getFlag());

			LCCollection* col_MCP = evtP->getCollection( _inputCollectionName ) ;
			int _nMCP=col_MCP->getNumberOfElements();

			std::vector<MCParticle* > FSParticles;
			std::vector<MCParticle* > FirstGenerationParticles; 

			int NParent   = 0; 
			int NDaughter = 0; 
			int mcPID     = 0;
			int Status    = -999;

			//TVector3 VTX, EndP; 

			for(int i0 = 0; i0 < _nMCP; i0++)
			{
				MCParticle *a_mcp = dynamic_cast<EVENT::MCParticle *>(col_MCP->getElementAt(i0));
				TVector3 curP = a_mcp->getMomentum();
				NParent = a_mcp->getParents().size();
				Status =  a_mcp->getGeneratorStatus();
				NDaughter = a_mcp->getDaughters().size();
				TVector3 mc3V = TVector3(a_mcp->getMomentum());
				//VTX = a_mcp->getVertex();
				//EndP = a_mcp->getEndpoint(); 
				double costheta = mc3V.CosTheta() ;  
				mcPID = a_mcp->getPDG();
				if ( Status!=1      ) continue; 
				if ( curP.Mag()<_momentumCut ) continue;
				if ( fabs(costheta)>0.98500  ) continue;	
				if(NParent == 0)
				{
					FirstGenerationParticles.push_back(a_mcp);
				}
				//if ( Status==1 && NDaughter == 0 ) 
				if ( Status==1 ) 
				{
					if( _rejectNeutrino>0 && (abs(mcPID)==12 || abs(mcPID)==14 || abs(mcPID)==16))  continue; // neutrinos cannot be detected
					FSParticles.push_back(a_mcp);
				}else{
					//printf("PDGID =  %10d, P = %103f\n", mcPID, curP.Mag());
				}
			}

			int NFSP = FSParticles.size();
			double ResoAmpli = 0; 	
			double CurrPAmpli = 0; 		
			double ScaledMomentum[3] = {0, 0, 0};
			TVector3 CurrPMomentum; 
			double ScaleFactor = 1.0; 
			double ValidChance = 1.0;
			double FlatRandom = 0; 

			for(int i1 = 0; i1 < NFSP; i1++)
			{
				ScaleFactor = 1.0; 
				MCParticle *b_mcp = FSParticles[i1];	//sort according to Energy?
				CurrPMomentum   = b_mcp->getMomentum();
				CurrPAmpli      = CurrPMomentum.Mag();
				//double costheta = CurrPMomentum.CosTheta();  
				_MCPEn  = b_mcp->getEnergy();
				_MCP[0] = b_mcp->getMomentum()[0];
				_MCP[1] = b_mcp->getMomentum()[1];
				_MCP[2] = b_mcp->getMomentum()[2];
				_MCPID  = b_mcp->getPDG();

				if(_EnableSmearing)
				{
					ResoAmpli = DetectorReso(b_mcp)[0];
				}

				if(_EnableEfficiency)
				{
					ValidChance = DetectorReso(b_mcp)[1];
					FlatRandom  = gRandom->Rndm(1);
					_FlatRnd    = FlatRandom;
				}

				if(FlatRandom < ValidChance)
				{
					if(_EnableSmearing)
					{
						if(abs(_MCPID) == 11)
						{
							TF1 *func = new TF1("func", CrystalBall, 0, 150, 4);
							func->SetParNames("#alpha", "n", "mean", "#sigma1");
							func->SetParameters(DetectorReso(b_mcp)[3], DetectorReso(b_mcp)[4], CurrPAmpli, ResoAmpli);
							ScaleFactor = 1.0*(func->GetRandom())/CurrPAmpli;
							delete func; 	
						}
						else if( abs(_MCPID) != 11 && fabs(b_mcp->getCharge())>0.01 ) 
						{
							ScaleFactor = 1.0*(CurrPAmpli+ResoAmpli)/CurrPAmpli;
						}
						else if(abs(_MCPID) == 22){
						  	ScaleFactor = 1.0*(CurrPAmpli+ResoAmpli)/CurrPAmpli;
						}

						//printf("PDGID =  %10d, P = %10.3f, factor = %10.3f\n", mcPID, CurrPAmpli, ScaleFactor);

						//Naive Protection
						if     (ScaleFactor <  0) ScaleFactor = 0.0;
						else if(ScaleFactor >  5) ScaleFactor = 5.0; 
						//if( fabs(costheta)>0.980) ScaleFactor *= ( 1.00 - (fabs(costheta)-0.980) / 0.020 ) ; 
					}
					_SF  = ScaleFactor; 
					_Eff = ValidChance; 
					ScaledMomentum[0] = CurrPMomentum.X()*ScaleFactor;
					ScaledMomentum[1] = CurrPMomentum.Y()*ScaleFactor;
					ScaledMomentum[2] = CurrPMomentum.Z()*ScaleFactor;

					ReconstructedParticleImpl * a_recoparticle = new ReconstructedParticleImpl();
					a_recoparticle->setEnergy(b_mcp->getEnergy()*ScaleFactor);
					a_recoparticle->setMomentum( ScaledMomentum );
					a_recoparticle->setCharge( b_mcp->getCharge() );
					a_recoparticle->setMass( b_mcp->getMass() );
					a_recoparticle->setType( b_mcp->getPDG() );				
					//
					TrackImpl                 * dummyTrack    = NULL;  
					ClusterImpl               * dummyCluster  = NULL;  
					if( fabs(b_mcp->getCharge())>0.01 ) {
						dummyTrack   = new TrackImpl() ; 
						a_recoparticle->addTrack  ( dummyTrack   ) ; // dummy track to make it look like a real particle !!! memory leakage
						_tracktrash.push_back( dummyTrack);
					}else{
						dummyCluster = new ClusterImpl(); 
						a_recoparticle->addCluster( dummyCluster ) ; // dummy cluster to make it look like a real  particle !!! memory leakage
						_clustertrash.push_back( dummyCluster);
					}
					//
					ReconstructedParticleImpl * dummyParticle = NULL;  
					dummyParticle   = new ReconstructedParticleImpl() ; 
					a_recoparticle->addParticle( dummyParticle ) ; // dummy cluster to make it look like a real  particle !!! memory leakage
					_particletrash.push_back( dummyParticle);
					//_particletrash.push_back( a_recoparticle);
					//
					fastrecoparticle->addElement(a_recoparticle);
					relNav.addRelation( a_recoparticle , b_mcp ) ;
					//
				}
			}
			fastrecoparticle->setDefault   ( true  ) ;
			fastrecoparticle->setSubset    ( false ) ;
			fastrecoparticle->setTransient ( true  ) ;
			evtP->addCollection( fastrecoparticle,             _recoParticleCollectionName );
			evtP->addCollection( relNav.createLCCollection() , _mcTruthCollectionName ) ;

			_nEvt ++ ;

		}		
		catch (lcio::DataNotAvailableException err) { }

	}  	  

}	

void PFAFastSimu::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
	FreeDelAll(_tracktrash);
	FreeDelAll(_clustertrash);
	FreeDelAll(_particletrash);
} 

void PFAFastSimu::end()
{
}
