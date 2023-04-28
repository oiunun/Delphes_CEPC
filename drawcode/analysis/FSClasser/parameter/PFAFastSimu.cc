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



#include <PFAFastSimu.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
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
const double PhotoReso = 0.16; 
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


}

float InvMass(float* P)
{
	float mass = 0; 

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

	//	printParameters();
	_Num = 0; 

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");	
	_outputTree->Branch("Num", &_Num, "Num/I");
	_outputTree->Branch("MCP", _MCP, "MCP[3]/F");
	_outputTree->Branch("MCPEn", &_MCPEn, "MCPEn/F");
	_outputTree->Branch("MCPID", &_MCPID, "MCPID/I");
	_outputTree->Branch("Eff", &_Eff, "Eff/F");
	_outputTree->Branch("FlatRnd", &_FlatRnd, "FlatRnd/F");
	_outputTree->Branch("SF", &_SF, "SF/F");	//Scale Factor
	_outputTree->Branch("ealpha", ealpha, "elapha[10]/D");	//Scale Factor

	ifstream electron("/scratchfs/bes/chenzx/cepc/FSPFA/para/electron.para");
	for(int i=0; i<10; i++)
	{
		electron>>ealpha[i]>>en[i]>>emean[i]>>esigma[i];
	}
	electron.close();

	//        ifstream muon("/scratchfs/bes/chenzx/cepc/FSPFA/para/muon.para");
	ifstream muon("/besfs/groups/higgs/users/chenzx/AnaGeo/fit/muon/1365/result.txt");
	//        ifstream muon("/besfs/groups/higgs/users/chenzx/AnaGeo/fit/muon/1465/result.txt");
	//        ifstream muon("/besfs/groups/higgs/users/chenzx/AnaGeo/fit/muon/1565/result.txt");
	//        ifstream muon("/besfs/groups/higgs/users/chenzx/AnaGeo/fit/muon/1665/result.txt");
	//        ifstream muon("/besfs/groups/higgs/users/chenzx/AnaGeo/fit/muon/1808/result.txt");
	for(int i=0; i<10; i++)
	{
		for(int j=0; j<19; j++)
		{
			muon>>mufrac[i][j]>>mumean1[i][j]>>musigma1[i][j]>>mumean2[i][j]>>musigma2[i][j]>>mueff[i][j];
		}
	}
	muon.close();


	ifstream pion("/scratchfs/bes/chenzx/cepc/FSPFA/para/pion.para");
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

			_eventNr=evtP->getEventNumber();

			//std::cout<<"Num "<<_Num<<std::endl;

			LCCollection *fastrecoparticle = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

			Cluflag.setBit(LCIO::CHBIT_LONG);

			fastrecoparticle->setFlag(Cluflag.getFlag());

			LCCollection* col_MCP = evtP->getCollection( "MCParticle" ) ;
			int _nMCP=col_MCP->getNumberOfElements();

			std::vector<MCParticle* > FSParticles;
			std::vector<MCParticle* > FirstGenerationParticles; 

			int NParent = 0; 
			int NDaughter = 0; 
			int mcPID = -999;

			for(int i0 = 0; i0 < _nMCP; i0++)
			{
				MCParticle *a_mcp = dynamic_cast<EVENT::MCParticle *>(col_MCP->getElementAt(i0));
				TVector3 curP = a_mcp->getMomentum();
				NParent = a_mcp->getParents().size();
				NDaughter = a_mcp->getDaughters().size();
				mcPID = a_mcp->getPDG();
				if(NParent == 0)
				{
					FirstGenerationParticles.push_back(a_mcp);
				}
				if(NDaughter == 0) 
				{
					if(abs(mcPID)==12 || abs(mcPID)==14 || abs(mcPID)==16)  continue; // neutrinos cannot be detected
					FSParticles.push_back(a_mcp);
				}
			}

			int NFSP = FSParticles.size();
			float ResoAmpli = 0; 	
			double CurrPAmpli = 0; 		
			float ScaleFactor = 1.0; 
			double ScaledMomentum[3] = {0, 0, 0};
			TVector3 CurrPMomentum; 
			float ValidChance = 1.0;
			double FlatRandom = 0; 

			for(int i1 = 0; i1 < NFSP; i1++)
			{

				MCParticle *b_mcp = FSParticles[i1];	//sort according to Energy?
				CurrPMomentum = b_mcp->getMomentum();
				CurrPAmpli = CurrPMomentum.Mag();
				_MCPEn = b_mcp->getEnergy();
				_MCP[0] = b_mcp->getMomentum()[0];
				_MCP[1] = b_mcp->getMomentum()[1];
				_MCP[2] = b_mcp->getMomentum()[2];
				_MCPID = b_mcp->getPDG();

				if(_EnableSmearing)
				{
					ResoAmpli = DetectorReso(b_mcp)[0];
				}

				if(_EnableEfficiency)
				{
					ValidChance = DetectorReso(b_mcp)[1];
					FlatRandom = gRandom->Rndm(1);
					_FlatRnd = FlatRandom;
				}

				if(FlatRandom < ValidChance)
				{
					if(abs(_MCPID) == 11)
					{
						TF1 *func = new TF1("func", CrystalBall, 0, 150, 4);
						func->SetParNames("#alpha", "n", "mean", "#sigma1");
						func->SetParameters(DetectorReso(b_mcp)[3], DetectorReso(b_mcp)[4], CurrPAmpli, ResoAmpli);
						ScaleFactor = 1.0*(func->GetRandom())/CurrPAmpli; 
					}
					else if(abs(_MCPID) == 13 || abs(_MCPID) == 211) 
					{
						ScaleFactor = 1.0*(CurrPAmpli+ResoAmpli)/CurrPAmpli;
					}
					else if(abs(_MCPID) == 22) ScaleFactor = 1.0*(CurrPAmpli+ResoAmpli)/CurrPAmpli;


					//Naive Protection
					if(ScaleFactor < 0) ScaleFactor = 1.0;
					else if(ScaleFactor > 10) ScaleFactor = 2.0; 

					_SF = ScaleFactor; 
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
					//					a_recoparticle->setPDG( b_mcp->getPDG() );				
					fastrecoparticle->addElement(a_recoparticle);
				}
			}
			evtP->addCollection( fastrecoparticle, "FSRecoParticle" );

			if(_Num % 2 == 0)
			{
				std::cout<<"Num "<<_Num<<std::endl; 
				_outputTree->Fill();
			}
			_Num ++;

		}		
		catch (lcio::DataNotAvailableException err) { }

	}  	  

}	

void PFAFastSimu::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;

	}

}



