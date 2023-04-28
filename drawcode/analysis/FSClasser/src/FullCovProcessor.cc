#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP

#include "marlin/FastMCParticleType.h"
#include "marlin/ErrorOfSigma.h"

//--- LCIO headers

#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "FullCovProcessor.h"
#include "cepcplotstyle.h"
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackerHitImpl.h>

using namespace lcio;

namespace marlin
{

	FullCovProcessor aFullCovProcessor;

	FullCovProcessor::FullCovProcessor() : Processor("FullCovProcessor"),
										   //_factory(NULL),
										   _nRun(-1),
										   _nEvt(-1)
	{

		// modify processor description
		_description = "FullCovProcessor creates ReconstrcutedParticles from MCParticles "
					   "according to the resolution given in the steering file.";

		// register steering parameters: name, description, class-variable, default value

		registerInputCollection(LCIO::MCPARTICLE,
								"InputCollectionName",
								"Name of the MCParticle input collection",
								_inputCollectionName,
								string("MCParticle"));

		registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
								 "RecoParticleCollectionName",
								 "Name of the ReconstructedParticles output collection",
								 _recoParticleCollectionName,
								 string("ReconstructedParticles"));

		registerOutputCollection(LCIO::LCRELATION,
								 "MCTruthMappingCollectionName",
								 "Name of the MCTruthMapping output collection",
								 _mcTruthCollectionName,
								 string("MCTruthMapping"));

		registerProcessorParameter("MomentumCut",
								   "No reconstructed particles are produced for smaller momenta (in [GeV])",
								   _momentumCut,
								   float(0.001));

		registerProcessorParameter("FTag",
								   "",
								   _Tag,
								   float(0));

		registerProcessorParameter("SaveNPZ",
								   "",
								   m_saveNPZ,
								   int(0));

		registerProcessorParameter("Scale",
								   " calibration parameter for the compensate of threshold ",
								   _Scale,
								   float(1.028));

		FloatVec chResDefault;
		chResDefault.push_back(0.005); //D0
		chResDefault.push_back(0.010);
		chResDefault.push_back(0.000);

		chResDefault.push_back(0.005); //Z0
		chResDefault.push_back(0.010);
		chResDefault.push_back(1.000);

		chResDefault.push_back(0.005); //T0

		chResDefault.push_back(0.001);	 //Theta
		chResDefault.push_back(0.001);	 //Phi
		chResDefault.push_back(0.00002);  //Momentum
		chResDefault.push_back(0.001);	 //Momentum

		chResDefault.push_back(3.00); // BField
		chResDefault.push_back(-1.0); //cosTheta_min
		chResDefault.push_back(1.0);  //cosTheta_max

		registerProcessorParameter("ChargedResolution",
								   "Resolution of charged particles in polar angle range:  d(1/P)  th_min  th_max",
								   _initChargedRes,
								   chResDefault,
								   chResDefault.size());

		FloatVec gammaResDefault;
		gammaResDefault.push_back(0.01);
		gammaResDefault.push_back(0.17);
		gammaResDefault.push_back(0.01);
		gammaResDefault.push_back(0.01);
		gammaResDefault.push_back(-0.98);
		gammaResDefault.push_back(0.98);

		registerProcessorParameter("PhotonResolution",
								   "Resolution dE/E=A+B/sqrt(E/GeV) of photons in polar angle range: A  B th_min  th_max",
								   _initPhotonRes,
								   gammaResDefault,
								   gammaResDefault.size());

		FloatVec hadronResDefault;
		hadronResDefault.push_back(0.04);
		hadronResDefault.push_back(0.50);
		hadronResDefault.push_back(0.01);
		hadronResDefault.push_back(0.01);
		hadronResDefault.push_back(-0.98);
		hadronResDefault.push_back( 0.98);

		registerProcessorParameter("NeutralHadronResolution",
								   "Resolution dE/E=A+B/sqrt(E/GeV) of neutral hadrons in polar angle range: A  B th_min  th_max",
								   _initNeutralHadronRes,
								   hadronResDefault,
								   hadronResDefault.size());

		registerProcessorParameter("MakePlots", "Make some plots for check", m_makeplots, 0);
		registerProcessorParameter("Smear", "modeling the detector res.", m_smear, 1);
		registerProcessorParameter("Perfect", "perfect detector, only MC ", m_perfect, 0);
		registerProcessorParameter("PID", "perfect PID ", m_pid, 0);
		registerProcessorParameter("Luxury", "save 4 momentum to ntuple ", m_luxury, 0);
		registerProcessorParameter("RejectNeutrino", "reject the undetectables  ", m_rejectNeutrino, 1);
	}

	void FullCovProcessor::init()
	{

		printParameters();

		_nRun = 0;
		_nEvt = 0;
		Nv = 4;
		Pm = 100;
		Ne = 0;
		dat.clear();

#ifdef MARLIN_CLHEP
		_ChResVec.clear();
		int index = 0;
		{
			double sgD0 = _initChargedRes[index++];
			double sgD0PtA = _initChargedRes[index++];
			double sgD0PtB = _initChargedRes[index++];

			double sgZ0 = _initChargedRes[index++];
			double sgZ0PtA = _initChargedRes[index++];
			double sgZ0PtB = _initChargedRes[index++];

			double sgT0 = _initChargedRes[index++];

			double sgPhi = _initChargedRes[index++];
			double sgTheta = _initChargedRes[index++];
			double sgPRelA = _initChargedRes[index++];
			double sgPRelB = _initChargedRes[index++];

			double bF = _initChargedRes[index++];
			double thMin = _initChargedRes[index++];
			double thMax = _initChargedRes[index++];

			_ChResVec.push_back(TrackResolution(
				sgD0,
				sgD0PtA,
				sgD0PtB,
				sgZ0,
				sgZ0PtA,
				sgZ0PtB,
				sgT0,
				sgPhi,
				sgTheta,
				sgPRelA,
				sgPRelB,
				bF, thMin, thMax));
         printf("In  reso = %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",              sgD0,             sgD0PtA,              sgD0PtB,              sgZ0,             sgZ0PtA,              sgZ0PtB,              sgPRelA,              sgPRelB );
         printf("Out reso = %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", _ChResVec[0].SgD0,_ChResVec[0].SgD0PtA, _ChResVec[0].SgD0PtB, _ChResVec[0].SgZ0,_ChResVec[0].SgZ0PtA, _ChResVec[0].SgZ0PtB, _ChResVec[0].SgPRelA, _ChResVec[0].SgPRelB );
		}
		_ClResVec.clear();
		{
			index = 0;
			double a = _initPhotonRes[index++];
			double b = _initPhotonRes[index++];
			double dthe = _initPhotonRes[index++];
			double dphi = _initPhotonRes[index++];
			double thMin = _initPhotonRes[index++];
			double thMax = _initPhotonRes[index++];
			_ClResVec.push_back(ClusterResolution(
				a,
				b,
				dthe,
				dphi,
				thMin, thMax));
		}
		{
			index = 0;
			double a = _initNeutralHadronRes[index++];
			double b = _initNeutralHadronRes[index++];
			double dthe = _initPhotonRes[index++];
			double dphi = _initPhotonRes[index++];
			double thMin = _initNeutralHadronRes[index++];
			double thMax = _initNeutralHadronRes[index++];
			_ClResVec.push_back(ClusterResolution(
				a,
				b,
				dthe,
				dphi,
				thMin, thMax));
		}

		streamlog_out(MESSAGE) << " FullCovProcessor::init() : registering FullCovParticleFactory " << endl;

#endif // MARLIN_CLHEP
		if (m_luxury > 0)
		{
			char *ntFullName = new char[100];
			sprintf(ntFullName, "nt%s", "MCTOPO");
			TTree *outputTree = new TTree(ntFullName, "");
			outputTree->SetAutoSave(32 * 1024 * 1024);
			m_ntp = new NTupleHelper(outputTree);
			delete ntFullName;
		}
		if (m_makeplots > 0)
		{
			SetPrelimStyle();
			h_Momentum[0] =  new TH1D("E_gamma", "Energy   of #gamma        ", 100, 0.0, 100.0);
			h_Momentum[1] =  new TH1D("p_electron", "Momentum of e^{+}      ", 180, 0.0, 180.0);
			h_Momentum[2] =  new TH1D("p_positron", "Momentum of e^{-}      ", 180, 0.0, 180.0);
			h_Momentum[3] =  new TH1D("p_muonplus", "Momentum of #mu^{+}    ", 180, 0.0, 180.0);
			h_Momentum[4] =  new TH1D("p_muonminus", "Momentum of #mu^{-}   ", 180, 0.0, 180.0);
			h_Momentum[17] = new TH1D("E_ISRgamma", "Energy   of ISRgamma   ", 100, 0.0, 100.0);

			h_Momentum[5] =  new TH1D("p_pionplus", "Momentum of #pi^{+}    ", 100, 0.0,  50.0);
			h_Momentum[6] =  new TH1D("p_pionminus", "Momentum of #pi^{-}   ", 100, 0.0,  50.0);
			h_Momentum[7] =  new TH1D("p_kaonplus", "Momentum of K^{+}      ", 100, 0.0,  50.0);
			h_Momentum[8] =  new TH1D("p_kaonminus", "Momentum of K^{-}     ", 100, 0.0,  50.0);
			h_Momentum[9] =  new TH1D("p_proton", "Momentum of p            ", 100, 0.0,  50.0);
			h_Momentum[10] = new TH1D("p_antiproton", "Momentum of #bar{p}  ", 100, 0.0,  50.0);
			h_Momentum[11] = new TH1D("p_neutron", "Momentum of n           ", 100, 0.0,  50.0);
			h_Momentum[12] = new TH1D("p_antineutron", "Momentum of #bar{n} ", 100, 0.0,  50.0);
			h_Momentum[13] = new TH1D("p_Klong", "Momentum of K_{L}         ", 100, 0.0,  50.0);
			h_Momentum[14] = new TH1D("p_pizero", "Momentum of #pi^{0}      ", 100, 0.0,  50.0);
			h_Momentum[15] = new TH1D("P_charged", "momenta  of charges     ", 100, 0.0, 100.0);
			h_Momentum[16] = new TH1D("P_neutral", "momenta  of neutrals    ", 100, 0.0, 100.0);
			h_Momentum[18] = new TH1D("K_Mom_in_B"   ,      "K  in B decays ", 100, 0.0,  50.0);
			h_Momentum[19] = new TH1D("K_Mom_in_D"   ,      "K  in D decays ", 100, 0.0,  50.0);
			h_Momentum[20] = new TH1D("Pi_Mom_in_B"  ,      "Pi in B decays ", 100, 0.0,  50.0);
			h_Momentum[21] = new TH1D("Pi_Mom_in_D"  ,      "Pi in D decays ", 100, 0.0,  50.0);
			//
			h_Mass[0] = new TH1D("M_pi0", "Mass of 2 #gamma    ", 100, 0.000, 0.400);
			h_Mass[1] = new TH1D("M_omega", "Mass of 2 #pi       ", 100, 0.700, 0.900);
			h_Mass[2] = new TH1D("M_Lambda", "Mass of p and #pi   ", 100, 1.050, 1.200);
			h_Mass[3] = new TH1D("M_phi", "Mass of KK          ", 100, 0.900, 1.200);
			h_Mass[4] = new TH1D("M_Jpsi", "Mass of di-lepton   ", 100, 2.100, 4.100);
			h_Mass[5] = new TH1D("M_Upsilon", "Mass of di-lepton   ", 100, 8.900, 9.900);
			h_Mass[6] = new TH1D("M_D0", "Mass of K #pi       ", 100, 1.600, 2.200);
			h_Mass[7] = new TH1D("M_D+", "Mass of K #pi #pi   ", 100, 1.600, 2.200);
			h_Mass[8] = new TH1D("M_Ds", "Mass of K K #pi     ", 100, 1.600, 2.200);
			h_Mass[9] = new TH1D("M_B0", "Mass of B0          ", 100, 4.200, 5.430);
			h_Mass[10] = new TH1D("M_B+", "Mass of B+          ", 100, 4.200, 5.430);
			h_Mass[11] = new TH1D("M_Bs", "Mass of Bs          ", 100, 4.800, 5.800);
			h_Mass[12] = new TH1D("M_Bc", "Mass of Bc          ", 100, 5.800, 6.800);
			//
			h_CosTheta[0] = new TH1D("all_particles", "Cos of Theta        ", 200, -1.0, 1.00);
			h_CosTheta[1] = new TH1D("ISR_Photon", "Cos of Theta        ", 200, -1.0, 1.00);
			//
			h_Ntrks[0] = new TH1D("Ntrks", "# of tracks         ", 150, 0.0, 150.00);
			h_Ntrks[1] = new TH1D("Ngamma", "# of photons        ", 150, 0.0, 150.00);
			h_Ntrks[2] = new TH1D("Nparticles", "# of particles      ", 150, 0.0, 150.00);
			h_Ntrks[3] = new TH1D("Nneutralhadron", "# of neutral hadons", 150, 0.0, 150.00);

			//
			h_Vertex[0] = new TH1D("D0", "Vertex of D0        ", 220, 0., 2200);
			h_Vertex[1] = new TH1D("D+", "Vertex of D+        ", 220, 0., 2200);
			h_Vertex[2] = new TH1D("Ds", "Vertex of Ds        ", 220, 0., 2200);
			h_Vertex[3] = new TH1D("B0", "Vertex of B0        ", 220, 0., 2200);
			h_Vertex[4] = new TH1D("B+", "Vertex of B+        ", 220, 0., 2200);
			h_Vertex[5] = new TH1D("Bs", "Vertex of Bs        ", 220, 0., 2200);
			h_Vertex[6] = new TH1D("Bc", "Vertex of Bc        ", 220, 0., 2200);

			h_ResD[0] = new TH1D("R0", "                    ", 200, -.04, .04);
			h_ResD[1] = new TH1D("Z0", "                    ", 200, -.1, .1);
			h_ResD[2] = new TH1D("Px", "                    ", 200, -.03, .03);
			h_ResD[3] = new TH1D("Py", "                    ", 200, -.03, .03);
			h_ResD[4] = new TH1D("Pp", "                    ", 200, -.003, .003);
			h_ResD[5] = new TH1D("En", "                    ", 200, -.01, .01);
		}
	}

	void FullCovProcessor::processRunHeader(LCRunHeader *run)
	{
		_nRun++;
		//FreeDelAll(_ptrash);
		//FreeDelAll(_tracktrash);
		//FreeDelAll(_clustertrash);
	}

	void FullCovProcessor::processEvent(LCEvent *evt)
	{
		_nEvt++;
		if (0 == _nEvt % 1000)
			printf("Run & event =  %8d\n", _nEvt);
		//

		//FreeDelAll(_ptrash);
		//FreeDelAll(_tracktrash);
		//FreeDelAll(_clustertrash);

		const LCCollection *mcpCol = evt->getCollection(_inputCollectionName);
		LCCollectionVec *recVec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		LCCollectionVec *parVec = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		LCCollectionVec *trkVec = new LCCollectionVec(LCIO::TRACK);
		LCCollectionVec *clsVec = new LCCollectionVec(LCIO::CLUSTER);

		//Cluflag.setBit(LCIO::CHBIT_LONG);
		//recVec->setFlag(Cluflag.getFlag());
		//parVec->setFlag(Cluflag.getFlag());
		//trkVec->setFlag(Cluflag.getFlag());
		//clsVec->setFlag(Cluflag.getFlag());

		LCRelationNavigator relNav(LCIO::RECONSTRUCTEDPARTICLE, LCIO::MCPARTICLE);

		vector<TLorentzVector> PhotonList, PiMinusList, PiPlusList, KPlusList, KMinusList, PList, AntiPList, ElectronList, PositronList, MuPlusList, MuMinusList;
		vector<MCParticle *> PhotonParent, PiMinusParent, PiPlusParent, KPlusParent, KMinusParent, PParent, AntiPParent, ElectronParent, PositronParent, MuPlusParent, MuMinusParent;

		vector<TLorentzVector> rawp4list;
		vector<double> pdgidlist;
		if (m_luxury > 0)
		{
			rawp4list.clear();
			pdgidlist.clear();
		}
		double VisEn = 0, ntrk = 0, ngam = 0, npar = 0, nNeuHad = 0;
		vector<double> npzEVT(Nv * Pm);
		if( m_saveNPZ > 0 ) {
			for (long unsigned int i = 0; i < Nv * Pm; i++)
			{
				npzEVT[i] = 0.0;
				dat.push_back(0.0);
			}
		}
		int bcl = 0;

		for (int i = 0; i < mcpCol->getNumberOfElements(); i++)
		{
			MCParticle *mcp = dynamic_cast<MCParticle *>(mcpCol->getElementAt(i));
			int status = (mcp)->getGeneratorStatus();
			int nParents = ((mcp)->getParents()).size();
			if (status == 1 || nParents > 0)
				continue;
			int pdg = mcp->getPDG();
			if (abs(pdg) < 6)
				bcl += abs(pdg);
		}
		bcl = (10 - bcl) / 2;
		if (bcl > 2)
			bcl = 2;

		if (bcl < 2)
			Ne++;
		int Np = 0;
		for ( long unsigned int i = 0; i < (long unsigned int) mcpCol->getNumberOfElements(); i++)
		{
			MCParticle *mcp = dynamic_cast<MCParticle *>(mcpCol->getElementAt(i));
			MCParticle *mot = NULL;

			int status = (mcp)->getGeneratorStatus();
			int nParents = ((mcp)->getParents()).size();
			int idParent =0, idGParent    = 0, nGGParents = 0;
			if ( nParents>0) {
				mot = mcp->getParents()[0];
				idGParent = ((mcp)->getParents()[0])->getPDG();
				nGGParents= ((mcp)->getParents()[0])->getParents().size();
			}	
			int idGGParent  = 0;
			if ( nGGParents>0) {
				idGGParent = (((mcp)->getParents()[0])->getParents()[0])->getPDG();
			}	

			if (m_rejectNeutrino != 0 &&
				(abs(mcp->getPDG()) == 12 || abs(mcp->getPDG()) == 14 || abs(mcp->getPDG()) == 16))
			{
				continue;
			}
			if (abs(mcp->getPDG()) == 22 && (nParents == 0 || (nParents > 0 && idParent == 22 && nGGParents == 0)))
			{
				continue;
			}

			TVector3 v(mcp->getMomentum());
			double En(mcp->getEnergy());
			double theta = v.Theta();
			double pmag = v.Mag();
			double mass = mcp->getMass();
			double charge = mcp->getCharge();
			//
			int pdgcode = mcp->getPDG();
			int pdgCODE = abs(pdgcode);
			double gabe = pmag / (mass + 1e-8);
			double vv = 0, xx[3] = {0, 0, 0};
			xx[0] = (mcp)->getEndpoint()[0];
			xx[1] = (mcp)->getEndpoint()[1];
			xx[2] = (mcp)->getEndpoint()[2];
			vv = pow(xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2], 0.5) / gabe * 1000;
			if (m_makeplots > 0)
			{
				if (pdgCODE == 421)
					h_Vertex[0]->Fill(vv, 1.0);
				if (pdgCODE == 411)
					h_Vertex[1]->Fill(vv, 1.0);
				if (pdgCODE == 431)
					h_Vertex[2]->Fill(vv, 1.0);
				if (pdgCODE == 511)
					h_Vertex[3]->Fill(vv, 1.0);
				if (pdgCODE == 521)
					h_Vertex[4]->Fill(vv, 1.0);
				if (pdgCODE == 531)
					h_Vertex[5]->Fill(vv, 1.0);
				if (pdgCODE == 541)
					h_Vertex[6]->Fill(vv, 1.0);
			}
			if (status == 1)
			{ // stable particles only

				if (!(m_perfect || pmag > _momentumCut))
					continue;
				if (m_makeplots > 0)
				{
					h_CosTheta[0]->Fill(cos(theta), 1.0);
					if (pdgcode == 22 && (nParents == 0 || (nParents > 0 && idParent == 22 && nGGParents == 0)))
						h_CosTheta[1]->Fill(cos(theta), 1.0);
					if (pdgcode == 211)
						h_Momentum[5]->Fill(pmag, 1.0);
					if (pdgcode == -211)
						h_Momentum[6]->Fill(pmag, 1.0);
					if (pdgcode == 321)
						h_Momentum[7]->Fill(pmag, 1.0);
					if (pdgcode == -321)
						h_Momentum[8]->Fill(pmag, 1.0);
					if (pdgcode == 2212)
						h_Momentum[9]->Fill(pmag, 1.0);
					if (pdgcode == -2212)
						h_Momentum[10]->Fill(pmag, 1.0);
					if (pdgcode == 2112)
						h_Momentum[11]->Fill(pmag, 1.0);
					if (pdgcode == -2112)
						h_Momentum[12]->Fill(pmag, 1.0);
					if (pdgcode == 130)
						h_Momentum[13]->Fill(pmag, 1.0);
					if (pdgcode == 111)
						h_Momentum[14]->Fill(pmag, 1.0);
				}
				if (m_makeplots > 0)
				{
					if (pdgcode == 22)
						h_Momentum[0]->Fill(pmag, 1.0);
					if (pdgcode == 11)
						h_Momentum[1]->Fill(pmag, 1.0);
					if (pdgcode == -11)
						h_Momentum[2]->Fill(pmag, 1.0);
					if (pdgcode == 13)
						h_Momentum[3]->Fill(pmag, 1.0);
					if (pdgcode == -13)
						h_Momentum[4]->Fill(pmag, 1.0);
					if (pdgcode != 22 && (nParents > 0 && charge != 0))
						h_Momentum[15]->Fill(pmag, 1.0);
					if (pdgcode != 22 && (nParents > 0 && charge == 0))
						h_Momentum[16]->Fill(pmag, 1.0);
					if (pdgcode == 22 && (nParents == 0 || (nParents > 0 && idParent == 22)))
						h_Momentum[17]->Fill(pmag, 1.0);

               Int_t bc1 = (abs(idGParent )/100)%10; 
               Int_t bc2 = (abs(idGGParent)/100)%10; 

               if(abs(pdgcode) ==  321 ) {
                  if ( (nParents>0 && bc1==5) ||(nGGParents>0 && bc2 ==5 ) ) 
                     h_Momentum[18]->Fill(pmag,1.0);
               }	
               if(abs(pdgcode) ==  321 ) {
                  if ( (nParents>0 && bc1==4) ||(nGGParents>0 && bc2 ==4 ) ) 
                     h_Momentum[19]->Fill(pmag,1.0);
               }	
               if(abs(pdgcode) ==  211 ) {
                  if ( (nParents>0 && bc1==5) ||(nGGParents>0 && bc2 ==5 ) ) 
                     h_Momentum[20]->Fill(pmag,1.0);
               }	
               if(abs(pdgcode) ==  211 ) {
                  if ( (nParents>0 && bc1==4) ||(nGGParents>0 && bc2 ==4 ) ) 
                     h_Momentum[21]->Fill(pmag,1.0);
               }	
            }


				ReconstructedParticleImpl *rec = new ReconstructedParticleImpl();
				TrackImpl *trk = new TrackImpl();
				ClusterImpl *cls = new ClusterImpl();
				TrackStateImpl *ts = new TrackStateImpl();

				int it = createTrack(trk, ts, mcp);
				int ic = createCluster(cls, mcp);
				int ip = createReconstructedParticle(rec, mcp, trk, cls);

				if (it == 0 && ic == 0 && ip == 0)
				{
					VisEn += En;
					npar++;
					double d0 = 0, z0 = 0;
					if (fabs(charge) > 0.01){
						d0 = rec->getTracks()[0]->getD0();
					}
					if (fabs(charge) > 0.01)
						ntrk++;
					if (pdgcode == 22)
						ngam++;
					if (pdgcode != 22 && fabs(charge) < 0.01)
						nNeuHad++;
					recVec->addElement(rec);
					relNav.addRelation(rec, mcp);
					if (fabs(charge) > 0.01)
						relNav.addRelation(rec->getTracks()[0], mcp);
					relNav.addRelation(rec->getClusters()[0], mcp);
					TLorentzVector p4(rec->getMomentum(), rec->getEnergy());
					//
					EVENT::TrackVec::const_iterator it_trk = (rec->getTracks()).begin();
					for (; it_trk != (rec->getTracks()).end(); it_trk++)
					{
						trkVec->addElement(*it_trk);
						_tracktrash.push_back(*it_trk);
					}
					EVENT::ClusterVec::const_iterator it_clu = (rec->getClusters()).begin();
					for (; it_clu != (rec->getClusters()).end(); it_clu++)
					{
						clsVec->addElement(*it_clu);
						_clustertrash.push_back(*it_clu);
					}
					//
					EVENT::ReconstructedParticleVec::const_iterator it_par = (rec->getParticles()).begin();
					for (; it_par != (rec->getParticles()).end(); it_par++)
					{
						parVec->addElement(*it_par);
						_ptrash.push_back(*it_par);
					}
					//
					if (m_luxury > 0)
					{
						rawp4list.push_back(p4);
						pdgidlist.push_back(pdgcode);
					}
					if( m_saveNPZ > 0 ) {
						if (i < Pm)
						{
							npzEVT[Np * Nv + 0] = p4.T();
							npzEVT[Np * Nv + 1] = p4.CosTheta();
							npzEVT[Np * Nv + 2] = p4.Phi();
							npzEVT[Np * Nv + 3] = pdgcode;
							//npzEVT[ Np*Nv + 4 ] =  d0      ;
						}
					}
					Np++;
					if (m_makeplots > 0)
					{
						if (pdgcode == 22)
						{
							PhotonList.push_back(p4);
							PhotonParent.push_back(mot);
						}
						if (pdgcode == 11)
						{
							ElectronList.push_back(p4);
							ElectronParent.push_back(mot);
						}
						if (pdgcode == -11)
						{
							PositronList.push_back(p4);
							PositronParent.push_back(mot);
						}
						if (pdgcode == -13)
						{
							MuPlusList.push_back(p4);
							MuPlusParent.push_back(mot);
						}
						if (pdgcode == 13)
						{
							MuMinusList.push_back(p4);
							MuMinusParent.push_back(mot);
						}
						if (pdgcode == 211)
						{
							PiPlusList.push_back(p4);
							PiPlusParent.push_back(mot);
						}
						if (pdgcode == -211)
						{
							PiMinusList.push_back(p4);
							PiMinusParent.push_back(mot);
						}
						if (pdgcode == 321)
						{
							KPlusList.push_back(p4);
							KPlusParent.push_back(mot);
						}
						if (pdgcode == -321)
						{
							KMinusList.push_back(p4);
							KMinusParent.push_back(mot);
						}
						if (pdgcode == 2212)
						{
							PList.push_back(p4);
							PParent.push_back(mot);
						}
						if (pdgcode == -2212)
						{
							AntiPList.push_back(p4);
							AntiPParent.push_back(mot);
						}
					}
				}
				//delete trk;
				//delete cls;
				//delete rec;
			}
		}
		//printf("VisEn = %10.2f\n",VisEn);
		recVec->setDefault(true);
		recVec->setSubset(false);
		recVec->setTransient(false);
		parVec->setDefault(false);
		parVec->setSubset(false);
		parVec->setTransient(false);
		trkVec->setDefault(false);
		trkVec->setSubset(false);
		trkVec->setTransient(false);
		clsVec->setDefault(false);
		clsVec->setSubset(false);
		clsVec->setTransient(false);

		evt->addCollection(recVec, _recoParticleCollectionName);
		evt->addCollection(parVec, "ParticleCollection");
		evt->addCollection(trkVec, "TrackCollection");
		evt->addCollection(clsVec, "ClusterCollection");
		evt->addCollection(relNav.createLCCollection(), _mcTruthCollectionName);

		//delete recVec;
		//delete parVec;
		//delete trkVec;
		//delete clsVec;
		//
		if (m_makeplots > 0)
		{
			//printf("No of Photon is %4d\n", PhotonList.size());
			//printf("No of Kaon+  is %4d\n", KPlusList.size());
			//printf("No of Kaon-  is %4d\n", KMinusList.size());
			h_Ntrks[0]->Fill(ntrk, 1.0);
			h_Ntrks[1]->Fill(ngam, 1.0);
			h_Ntrks[2]->Fill(npar, 1.0);
			h_Ntrks[3]->Fill(nNeuHad, 1.0);
			// pi0
			if (PhotonList.size() > 1)
			{
				for (unsigned int i = 0; i < PhotonList.size() - 1; i++)
				{
					if ( PhotonParent.size()> i && PhotonParent[i]->getPDG() != 111)
						continue;
					for (unsigned int j = i + 1; j < PhotonList.size(); j++)
					{
						if (PhotonParent[j] != PhotonParent[i])
							continue;
						h_Mass[0]->Fill((PhotonList[i] + PhotonList[j]).M());
					}
				}
			}
			// omega
			if (PiPlusList.size() > 0 && PiMinusList.size() > 0)
			{
				for (unsigned int i = 0; i < PiPlusList.size(); i++)
				{
					if (PiPlusParent[i]->getPDG() != 223)
						continue;
					for (unsigned int j = 0; j < PiMinusList.size(); j++)
					{
						if (PiMinusParent[j] != PiPlusParent[i])
							continue;
						h_Mass[1]->Fill((PiPlusList[i] + PiMinusList[j]).M());
					}
				}
			}
			// Lambda
			if (PiPlusList.size() > 0 && AntiPList.size() > 0)
			{
				for (unsigned int i = 0; i < PiPlusList.size(); i++)
				{
					if (PiPlusParent[i]->getPDG() != -3122)
						continue;
					for (unsigned int j = 0; j < AntiPList.size(); j++)
					{
						if (AntiPParent[j] != PiPlusParent[i])
							continue;
						h_Mass[2]->Fill((PiPlusList[i] + AntiPList[j]).M());
					}
				}
			}

			if (PiMinusList.size() > 0 && PList.size() > 0)
			{
				for (unsigned int i = 0; i < PiMinusList.size(); i++)
				{
					if (PiMinusParent[i]->getPDG() != 3122)
						continue;
					for (unsigned int j = 0; j < PList.size(); j++)
					{
						if (PParent[j] != PiMinusParent[i])
							continue;
						h_Mass[2]->Fill((PiMinusList[i] + PList[j]).M());
					}
				}
			}
			// phi
			if (KPlusList.size() > 0 && KMinusList.size() > 0)
			{
				for (unsigned int i = 0; i < KPlusList.size(); i++)
				{
					int ip1 = KPlusParent[i]->getPDG();
					for (unsigned int j = 0; j < KMinusList.size(); j++)
					{
						int ip2 = KMinusParent[j]->getPDG();
						if (ip1 == 333 && KPlusParent[i] == KMinusParent[j])
							h_Mass[3]->Fill((KPlusList[i] + KMinusList[j]).M());
					}
				}
			}
			//J/psi and Upsilon
			if (MuPlusList.size() > 0 && MuMinusList.size() > 0)
			{
				for (unsigned int i = 0; i < MuPlusList.size(); i++)
				{
					int ip1 = MuPlusParent[i]->getPDG();
					for (unsigned int j = 0; j < MuMinusList.size(); j++)
					{
						int ip2 = MuMinusParent[j]->getPDG();
						if (ip1 == 443 && ip2 == 443 && MuPlusParent[i] == MuMinusParent[j])
							h_Mass[4]->Fill((MuPlusList[i] + MuMinusList[j]).M());
						if (ip1 == 553 && ip2 == 553 && MuPlusParent[i] == MuMinusParent[j])
							h_Mass[5]->Fill((MuPlusList[i] + MuMinusList[j]).M());
					}
				}
			}
			if (ElectronList.size() > 0 && PositronList.size() > 0)
			{
				for (unsigned int i = 0; i < ElectronList.size(); i++)
				{
					int ip1 = ElectronParent[i]->getPDG();
					for (unsigned int j = 0; j < PositronList.size(); j++)
					{
						int ip2 = PositronParent[j]->getPDG();
						if (ip1 == 443 && ip2 == 443 && ElectronParent[i] == PositronParent[j])
							h_Mass[4]->Fill((ElectronList[i] + PositronList[j]).M());
						if (ip1 == 553 && ip2 == 553 && ElectronParent[i] == PositronParent[j])
							h_Mass[5]->Fill((ElectronList[i] + PositronList[j]).M());
					}
				}
			}
			// D0
			if (KPlusList.size() > 0 && PiMinusList.size() > 0)
			{
				for (unsigned int i = 0; i < KPlusList.size(); i++)
				{
					int ip1 = abs(KPlusParent[i]->getPDG());
					for (unsigned int j = 0; j < PiMinusList.size(); j++)
					{
						int ip2 = abs(PiMinusParent[j]->getPDG());
						if (ip1 == 421 && ip2 == 421 && KPlusParent[i] == PiMinusParent[j])
							h_Mass[6]->Fill((KPlusList[i] + PiMinusList[j]).M());
					}
				}
			}
			if (KMinusList.size() > 0 && PiPlusList.size() > 0)
			{
				for (unsigned int i = 0; i < KMinusList.size(); i++)
				{
					int ip1 = abs(KMinusParent[i]->getPDG());
					for (unsigned int j = 0; j < PiPlusList.size(); j++)
					{
						int ip2 = abs(PiPlusParent[j]->getPDG());
						if (ip1 == 421 && ip2 == 421 && KMinusParent[i] == PiPlusParent[j])
							h_Mass[6]->Fill((KMinusList[i] + PiPlusList[j]).M());
					}
				}
			}

			// D+/-
			if (KPlusList.size() > 0 && PiMinusList.size() > 1)
			{
				for (unsigned int i = 0; i < KPlusList.size(); i++)
				{
					int ip1 = abs(KPlusParent[i]->getPDG());
					for (unsigned int j = 0; j < PiMinusList.size(); j++)
					{
						int ip2 = abs(PiMinusParent[j]->getPDG());
						for (unsigned int k = j + 1; k < PiMinusList.size(); k++)
						{
							int ip3 = abs(PiMinusParent[k]->getPDG());
							if (ip1 == 411 && ip2 == 411 && ip3 == 411 && KPlusParent[i] == PiMinusParent[j] && KPlusParent[i] == PiMinusParent[k])
								h_Mass[7]->Fill((KPlusList[i] + PiMinusList[j] + PiMinusList[k]).M());
						}
					}
				}
			}
			if (KMinusList.size() > 0 && PiPlusList.size() > 1)
			{
				for (unsigned int i = 0; i < KMinusList.size(); i++)
				{
					int ip1 = abs(KMinusParent[i]->getPDG());
					for (unsigned int j = 0; j < PiPlusList.size(); j++)
					{
						int ip2 = abs(PiPlusParent[j]->getPDG());
						for (unsigned int k = j + 1; k < PiPlusList.size(); k++)
						{
							int ip3 = abs(PiPlusParent[k]->getPDG());
							if (ip1 == 411 && ip2 == 411 && ip3 == 411 && KMinusParent[i] == PiPlusParent[j] && KMinusParent[i] == PiPlusParent[k])
								h_Mass[7]->Fill((KMinusList[i] + PiPlusList[j] + PiPlusList[k]).M());
						}
					}
				}
			}
			// Ds
			if (KPlusList.size() > 0 && KMinusList.size() > 0 && PiPlusList.size() > 0)
			{
				for (unsigned int i = 0; i < KPlusList.size(); i++)
				{
					int ip1 = abs(KPlusParent[i]->getPDG());
					for (unsigned int j = 0; j < KMinusList.size(); j++)
					{
						int ip2 = abs(KMinusParent[j]->getPDG());
						for (unsigned int k = j + 1; k < PiPlusList.size(); k++)
						{
							int ip3 = abs(PiPlusParent[k]->getPDG());
							if (ip1 == 431 && ip2 == 431 && ip3 == 431 && KPlusParent[i] == KMinusParent[j] && KPlusParent[i] == PiPlusParent[k])
								h_Mass[8]->Fill((KPlusList[i] + KMinusList[j] + PiPlusList[k]).M());
						}
					}
				}
			}
			if (KPlusList.size() > 0 && KMinusList.size() > 0 && PiMinusList.size() > 0)
			{
				for (unsigned int i = 0; i < KPlusList.size(); i++)
				{
					int ip1 = abs(KPlusParent[i]->getPDG());
					for (unsigned int j = 0; j < KMinusList.size(); j++)
					{
						int ip2 = abs(KMinusParent[j]->getPDG());
						for (unsigned int k = j + 1; k < PiMinusList.size(); k++)
						{
							int ip3 = abs(PiMinusParent[k]->getPDG());
							if (ip1 == 431 && ip2 == 431 && ip3 == 431 && KPlusParent[i] == KMinusParent[j] && KPlusParent[i] == PiMinusParent[k])
								h_Mass[8]->Fill((KPlusList[i] + KMinusList[j] + PiMinusList[k]).M());
						}
					}
				}
			}
		}
		//
		if (m_luxury > 0)
		{
			m_ntp->fill4Momentum("idx1", "raw_", rawp4list, rawp4list.size());
			m_ntp->fillArray("idx2", "pdgid", pdgidlist, pdgidlist.size());
			m_ntp->write();
		}
		if (bcl < 2)
		{
			if( m_saveNPZ > 0 ) {
				tag.push_back(bcl);
				for ( long unsigned int i = 0; i < Nv * Pm; i++)
				{
					dat[Nv * Pm * (Ne - 1) + i] = npzEVT[i];
				}
			}
		}
		npzEVT.clear();
	}

	void FullCovProcessor::check(LCEvent *evt)
	{
	}

	int FullCovProcessor::createReconstructedParticle(ReconstructedParticleImpl *rec, const MCParticle *mcp, Track *trk, Cluster *cls)
	{
		float m = mcp->getMass();
		rec->setType(mcp->getPDG());
		float p[3], e;
		if (fabs(mcp->getCharge()) > 0.01)
		{
			const TrackState *ts = trk->getTrackState(TrackState::AtIP);
			if (ts)
			{
				double p0[3] = {mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2]};
				double pm0(pow(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2], 0.5));
				float pt(_ChResVec[0].BF * 2.99792e-4 / fabs(ts->getOmega()));
				p[0] = pt * cos(ts->getPhi());
				p[1] = pt * sin(ts->getPhi());
				p[2] = pt / ts->getTanLambda();
				double pm(pow(p[0] * p[0] + p[1] * p[1] + p[2] * p[2], 0.5));
				if ( m_pid == 0 )
					m = 0.13957;
				e = pow(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + m * m, 0.5);
				//
				if (m_makeplots > 0)
				{
					h_ResD[2]->Fill((p[0] - p0[0]) / pm0, 1.0);
					h_ResD[3]->Fill((p[1] - p0[1]) / pm0, 1.0);
					h_ResD[4]->Fill((pm - pm0) / pm0, 1.0);
					h_ResD[5]->Fill((e - mcp->getEnergy()) / mcp->getEnergy(), 1.0);
					//printf("Energy =  %10.3f %10.3f %10.3f\n", e, pow(pt*pt+p[2]*p[2], 0.5), m );
				}
				rec->addTrack(trk);
			}
			else
			{
				return 1;
			}
		}
		else
		{
			if ( m_pid == 0 ) m = 0;
			const Hep3Vector v(mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2]);
			const double theta = v.theta();
			const double phi = v.phi();
			pair<double, double> angles = ensureThetaBounds(phi, theta);
			Hep3Vector dir;
			e = cls->getEnergy();
			if (mcp->getPDG() == 22)
			{
				const double deltaPhi = RandGauss::shoot(0.0, _ClResVec[0].dPhi);
				const double deltaThe = RandGauss::shoot(0.0, _ClResVec[0].dThe);
				if (!m_perfect && m_smear)
					angles = ensureThetaBounds(phi + deltaPhi, theta + deltaThe);
				dir = Hep3Vector(
					sin(angles.second) * cos(angles.first),
					sin(angles.second) * sin(angles.first),
					cos(angles.second));
			}
			else
			{
				const double deltaPhi = RandGauss::shoot(0.0, _ClResVec[1].dPhi);
				const double deltaThe = RandGauss::shoot(0.0, _ClResVec[1].dThe);
				if (!m_perfect && m_smear)
					angles = ensureThetaBounds(phi + deltaPhi, theta + deltaThe);
				dir = Hep3Vector(
					sin(angles.second) * cos(angles.first),
					sin(angles.second) * sin(angles.first),
					cos(angles.second));
			}
			dir.setMag(e);
			p[0] = dir.x();
			p[1] = dir.y();
			p[2] = dir.z();
		}
		rec->setMomentum(p);
		rec->setEnergy(e);
		rec->setMass(m);
		rec->setCharge(mcp->getCharge());
		rec->addCluster(cls);
		ReconstructedParticleImpl *dummy = new ReconstructedParticleImpl();
		dummy->id();
		rec->addParticle(dummy); // dummy track to make it look like a real particle !!! memory leakage
		//delete dummy;
		return 0;
	}

	int FullCovProcessor::createTrack(TrackImpl *trk, TrackStateImpl *ts, const MCParticle *mcp)
	{
		const Hep3Vector v(mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2]);
		const double costheta = cos(v.theta());
		const double sintheta = sin(v.theta());
		pair<double, double> angles = ensureThetaBounds(v.phi(), v.theta());
		//printf("reso = %8.3f %8.3f %8.3f \n", _ChResVec[0].SgD0, _ChResVec[0].SgD0PtA, _ChResVec[0].SgD0PtB);
		if (m_perfect || (costheta < _ChResVec[0].ThMax && costheta > _ChResVec[0].ThMin))
		{
			trk->id();
			const double pt = v.perp();
			const double p = v.mag();
			const double theta = v.theta();
			const double phi = v.phi();
			//
			const double SigD0 = _ChResVec[0].SgD0 + _ChResVec[0].SgD0PtA * exp(-1.0 * abs(_ChResVec[0].SgD0PtB) * pt);
			const double SigZ0 = _ChResVec[0].SgZ0 + _ChResVec[0].SgZ0PtA * exp(-1.0 * abs(_ChResVec[0].SgZ0PtB) * pt);
			// shortcuts for other resolutions
			const double SigT0 = _ChResVec[0].SgT0;
			const double SigPhi = _ChResVec[0].SgPhi;
			const double SigTheta = _ChResVec[0].SgTheta;
			const double SigU = SigD0;
			const double SigV = SigZ0 * sintheta;
			// draw random noise
			const float deltaD0 = SigD0 * RandGauss::shoot(0.0, 1.0);
			const float deltaZ0 = SigZ0 * RandGauss::shoot(0.0, 1.0);
			const float deltaT0 = SigT0 * RandGauss::shoot(0.0, 1.0);
			const float deltaPhi = SigPhi * RandGauss::shoot(0.0, 1.0);
			const float deltaTheta = SigTheta * RandGauss::shoot(0.0, 1.0);
			// smear direction angles phi,theta ensuring correct bounds
			if (!m_perfect && m_smear)
				angles = ensureThetaBounds(phi + deltaPhi, theta + deltaTheta);
			//
			const double relPt = pow(pow(_ChResVec[0].SgPRelA * p, 2) + pow(_ChResVec[0].SgPRelB, 2), 0.5) / _ChResVec[0].BF;
			//if( p> 10) printf("2 p relative = %8.2f %10.5f %10.5f\n", p, _ChResVec[0].SgPRelA,   _ChResVec[0].SgPRelB/2.0/pt/pow(sintheta,0.5));
			//if( p> 10) printf("3 p relative = %8.2f %10.5f %10.5f\n", p, _ChResVec[0].SgPRelA,   _ChResVec[0].SgPRelB/3.0/pt/pow(sintheta,0.5));
			const double SigPt = relPt * pt;
			const double SigP = SigPt / sin(angles.second);
			// var(q/p) = (d(1/p)/dp)² * var(p) = (-1/p²)² * var(p)
			const double SigQOP = SigPt / (pt * pt);
			// smear momentum
			const float deltaP = SigP * RandGauss::shoot(0.0, 1.0);

			// smear the position
			Hep3Vector pos0(mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2]);
			Hep3Vector pos(deltaD0 * sin(phi), -deltaD0 * cos(phi), deltaZ0);
			if (m_makeplots > 0)
			{
				h_ResD[0]->Fill(pos[0], 1.0);
				h_ResD[0]->Fill(pos[1], 1.0);
				h_ResD[1]->Fill(pos[2], 1.0);
			}
			if (!m_perfect && m_smear)
				pos = pos + pos0;
			float D0 = pow(pos[0] * pos[0] + pos[1] * pos[1], 0.5);
			float Z0 = pos[2];
			//
			// compute smeared direction vector
			const Hep3Vector dir(
				sin(angles.second) * cos(angles.first),
				sin(angles.second) * sin(angles.first),
				cos(angles.second));
			// compute smeared momentum vector
			Hep3Vector mom = p * dir;
			if (!m_perfect && m_smear)
				mom = (p + deltaP) * dir;
			if (m_makeplots > 0)
			{
			}
			//if( p> 30) printf("p relative = %10.4f %10.4f %10.4f\n", SigP, SigP/p*100, p+deltaP);

			EVENT::FloatVec covMatrix;

			covMatrix.resize(15);

			for (unsigned icov = 0; icov < covMatrix.size(); ++icov)
			{
				covMatrix[icov] = 0;
			}

			covMatrix[0] = (SigD0 * SigD0);							  //sigma_d0^2
			covMatrix[2] = (SigPhi * SigPhi);						  //sigma_phi0^2
			covMatrix[5] = (SigQOP * SigQOP);						  //sigma_omega^2
			covMatrix[9] = (SigZ0 * SigZ0);							  //sigma_z0^2
			covMatrix[14] = (SigTheta * SigTheta) / pow(costheta, 4); //sigma_tanl^2

			float Omega = (_ChResVec[0].BF * 2.99792e-4 * mcp->getCharge()) / ((p + deltaP) * sin(angles.second));
			//float  Omega  =   ((p+deltaP)*sin(angles.second)) * mcp->getCharge() * _ChResVec[0].BF * 2.99792e-4 ;
			float TanLambda = tan(angles.second);
			float vtx[3];
			vtx[0] = mcp->getVertex()[0];
			vtx[1] = mcp->getVertex()[1];
			vtx[2] = mcp->getVertex()[2];

			ts->setD0(D0);
			ts->setZ0(Z0);
			ts->setPhi(angles.first);
			ts->setOmega(Omega);
			ts->setTanLambda(TanLambda);
			ts->setCovMatrix(covMatrix);
			ts->setReferencePoint(vtx);
			ts->setLocation(TrackState::AtIP);
			//
			trk->addTrackState(ts);
			trk->subdetectorHitNumbers().resize(10);
		}
		else
		{
			return 1;
		}
		return 0;
	}

	int FullCovProcessor::createCluster(ClusterImpl *cls, const MCParticle *mcp)
	{

		const HepLorentzVector v(mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], mcp->getEnergy());

		// find resolution for polar angle
		double theta = v.theta();
		double phi = v.theta();
		double costheta = cos(theta);

		pair<double, double> resolution = make_pair(-1., -1.);

		if (m_perfect || (costheta <= _ClResVec[0].ThMax && costheta > _ClResVec[0].ThMin))
		{
			if (mcp->getPDG() == 22)
			{
				resolution.first = _ClResVec[0].A;
				resolution.second = _ClResVec[0].B;
			}
		}
		if (m_perfect || (costheta <= _ClResVec[1].ThMax && costheta > _ClResVec[1].ThMin))
		{
			if (abs(mcp->getPDG()) != 22)
			{
				resolution.first = _ClResVec[1].A;
				resolution.second = _ClResVec[1].B;
			}
		}

		double e = 0;

		if (resolution.first > -1e-10)
		{

			double E = v.e();

			if (!m_perfect && m_smear)
			{
				double Eres = sqrt(resolution.first * resolution.first + resolution.second * resolution.second / E);
				E = E * RandGauss::shoot(1.0, Eres);
			}
			e = fabs(E);
		}
		cls->setEnergy(e);
		return 0;
	}

	void FullCovProcessor::end()
	{

		streamlog_out(MESSAGE4) << "FullCovProcessor::end()  " << name()
								<< " processed " << _nEvt << " events in " << _nRun << " runs "
								<< endl;
		//cout<<" # of entries  " <<dat.size()<<" "<<tag.size() <<endl;

		if( m_saveNPZ > 0 ) {
			cnpy::npz_save("out.npz","X",&dat[0],{Ne,  Pm,  Nv}, "w");
			cnpy::npz_save("out.npz","Y",&tag[0],{Ne          }, "a");
		}
		dat.clear();
		tag.clear();
		if (m_luxury > 0)
		{
			delete m_ntp;
		}
		if (m_makeplots > 0)
		{
			for (int i = 0; i < 6; i++)
			{
				h_ResD[i]->Write();
				//
				TH1D *h1 = new TH1D("h1", "", h_ResD[i]->GetNbinsX(), h_ResD[i]->GetXaxis()->GetXmin(), h_ResD[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char *)"#Delta", (char *)"Entries");
				TH1D *h2 = 0;
				char filename[256], title[256], Error[256];
				//
				sprintf(filename, "figs/%s", (h_ResD[i]->GetName()));
				sprintf(title, "%s", (h_ResD[i]->GetTitle()));
				sprintf(Error, "#sigma_{m} = %8.3f #pm %6.3f\n", h_ResD[i]->GetStdDev() * 1000, h_ResD[i]->GetStdDevError() * 1000);
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_ResD[i], (char *)h_ResD[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, false, false, Error);
				//
				sprintf(filename, "figs/%s_log", (h_ResD[i]->GetName()));
				sprintf(title, "%s", (h_ResD[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_ResD[i], (char *)h_ResD[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, true, false, title);
				delete h1;
				delete h_ResD[i];
			}

			for (int i = 0; i < 22; i++)
			{
				h_Momentum[i]->Write();
				//
				TH1D *h1 = new TH1D("h1", "", h_Momentum[i]->GetNbinsX(), h_Momentum[i]->GetXaxis()->GetXmin(), h_Momentum[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char *)"P(GeV/c)", (char *)"Entries/1.0GeV/c");
				TH1D *h2 = 0;
				char filename[256], title[256];
				//
				sprintf(filename, "figs/%s", (h_Momentum[i]->GetName()));
				sprintf(title, "%s", (h_Momentum[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_Momentum[i], (char *)h_Momentum[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, false, false, title);
				//
				sprintf(filename, "figs/%s_log", (h_Momentum[i]->GetName()));
				sprintf(title, "%s", (h_Momentum[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_Momentum[i], (char *)h_Momentum[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, true, false, title);
				delete h1;
				delete h_Momentum[i];
			}
			char filename[256], title[256];
			char *names[] = {(char *)"", (char *)"MC Truth", (char *)"Exp Fit"};
			//
			for (int i = 0; i < 7; i++)
			{
				if (h_Vertex[i]->Integral() < 10)
					continue;
				NameAxes(h_Vertex[i], (char *)"Vertex(#mu m)", (char *)"Entries / 10 #mu m");
				h_Vertex[i]->Fit("expo", "+");
				sprintf(filename, "figs/Vertex_%s", h_Vertex[i]->GetName());
				sprintf(title, "%s", (h_Vertex[i]->GetTitle()));
				PlotDataFit(filename, h_Vertex[i], (char *)h_Vertex[i]->GetTitle(), names, true, 0.65, 0.60, 0.90, 0.80, title);
			}
			for (int i = 0; i < 7; i++)
			{
				h_Vertex[i]->Write();
				delete h_Vertex[i];
			}

			for (int i = 0; i < 13; i++)
			{
				//
				if (h_Mass[i]->Integral() < 10)
					continue;
				TH1D *h1 = new TH1D("h1", "", h_Mass[i]->GetNbinsX(), h_Mass[i]->GetXaxis()->GetXmin(), h_Mass[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char *)"Mass(GeV/c^{2})", (char *)"Entries");
				TH1D *h2 = NULL;
				char filename[256], title[256], Error[256];
				//
				sprintf(title, "%s", (h_Mass[i]->GetTitle()));
				sprintf(filename, "figs/%s", (h_Mass[i]->GetName()));
				sprintf(Error, "#sigma_{m} = %6.1f #pm %4.1f(MeV)\n", h_Mass[i]->GetStdDev() * 1000, h_Mass[i]->GetStdDevError() * 1000);
				PlotDataMC(filename,
						   h1, (char *)h_Mass[i]->GetTitle(),
						   h_Mass[i], (char *)h_Mass[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, false, false, Error);
				/*
				sprintf(filename,"figs/%s_log", (h_Mass[i]->GetName()));
				sprintf(title,"%s", (h_Mass[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Mass[i], (char*)h_Mass[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
					  );
				if( h_Mass[i]->Integral()>10){
					sprintf(title,"%s", (h_Mass[i]->GetTitle()));
					sprintf(filename,"figs/fit_%s", (h_Mass[i]->GetName()));
					h_Mass[i]->Fit("gaus");
					PlotDataFit(filename, h_Mass[i], (char*)h_Mass[i]->GetTitle(), names, true, 0.65,0.60,0.90,0.80, title);
				}
			   */
				delete h1;
			}
			for (int i = 0; i < 13; i++)
			{
				h_Mass[i]->Write();
				delete h_Mass[i];
			}
			for (int i = 0; i < 4; i++)
			{
				h_Ntrks[i]->Write();
				//
				TH1D *h1 = new TH1D("h1", "", h_Ntrks[i]->GetNbinsX(), h_Ntrks[i]->GetXaxis()->GetXmin(), h_Ntrks[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char *)"N_{particles}", (char *)"Entries/1.00");
				TH1D *h2 = 0;
				char filename[256], title[256];
				//
				sprintf(filename, "figs/%s", (h_Ntrks[i]->GetName()));
				sprintf(title, "%s", (h_Ntrks[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_Ntrks[i], (char *)h_Ntrks[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, false, false, title);
				//
				sprintf(filename, "figs/%s_log", (h_Ntrks[i]->GetName()));
				sprintf(title, "%s", (h_Ntrks[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_Ntrks[i], (char *)h_Ntrks[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, true, false, title);
				delete h1;
				delete h_Ntrks[i];
			}

			for (int i = 0; i < 2; i++)
			{
				h_CosTheta[i]->Write();
				//
				TH1D *h1 = new TH1D("h1", "", h_CosTheta[i]->GetNbinsX(), h_CosTheta[i]->GetXaxis()->GetXmin(), h_CosTheta[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char *)"cos#theta", (char *)"Entries/0.01");
				TH1D *h2 = 0;
				char filename[256], title[256];
				//
				sprintf(filename, "figs/%s", (h_CosTheta[i]->GetName()));
				sprintf(title, "%s", (h_CosTheta[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_CosTheta[i], (char *)h_CosTheta[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, false, false, title);
				//
				sprintf(filename, "figs/%s_log", (h_CosTheta[i]->GetName()));
				sprintf(title, "%s", (h_CosTheta[i]->GetTitle()));
				PlotDataMC(filename,
						   h1, (char *)"",
						   h_CosTheta[i], (char *)h_CosTheta[i]->GetTitle(),
						   h2, (char *)"",
						   h2, (char *)"",
						   h2, (char *)"",
						   true, true, false, title);
				delete h1;
				delete h_CosTheta[i];
			}
		}
	}
} // namespace marlin
//**********************************************
//**********************************************
int FullCovProcessor::FindParton(MCParticle *mcp)
{
	int ret = 0;

	if (mcp->getParents().size() == 0) // original particle, no parent
	{
		if (abs(mcp->getPDG()) <= 6)
		{
			return mcp->getPDG();
		}
		else
			return 0;
	}
	else if (mcp->getParents().size() == 1) // Only one parent
	{
		MCParticle *mcpa = mcp->getParents()[0];
		int ngp = mcpa->getParents().size();
		//
		if (ngp == 1)
		{
			if (mcp->getPDG() == mcpa->getPDG())
			{
				// mother is the same particle
				ret = FindParton(mcpa);
			}
			else if (abs(mcp->getPDG()) <= 6 && (abs(mcpa->getPDG()) >= 22 && abs(mcpa->getPDG()) <= 25))
			{
				// parent is gamma/gluon/W/Z and itself is a quark
				ret = mcp->getPDG();
			}
			else if (mcp->getPDG() == 21 && mcpa->getPDG() == 25)
			{
				// parent is Higgs and itself is a gluon
				ret = mcp->getPDG();
			}
			else if (abs(mcp->getPDG()) <= 6 && mcpa->getPDG() == 21)
			{
				// parent is gluon and itself is a quark, this kind of object marked by +100
				ret = 100 + mcp->getPDG();
			}
			else if (abs(mcpa->getPDG()) <= 6 && mcp->getPDG() == 21)
			{
				// parent is quark and itself is a gluon
				ret = FindParton(mcpa);
			}
		}
		else
		{
			if (abs(mcpa->getPDG()) >= 81 && abs(mcpa->getPDG()) <= 100)
			{
				// parent is intermediate particle; need to check the MANY(>1) grandparents
				if (abs(mcp->getPDG()) <= 6 || (abs(mcp->getPDG()) >= 11 && abs(mcp->getPDG()) <= 16))
				{
					// itself is a quark/lepton, jump over the intermediate paticle, search for the same quark/lepton, then ...
					for (unsigned int i = 0; i < mcpa->getParents().size(); i++)
					{
						MCParticle *mcpb = mcpa->getParents()[i];
						if (mcpb->getPDG() == mcp->getPDG())
						{
							//printf( "here: %5d vs. %5d (%5d)\n", mcp->getPDG(),mcpb->getPDG(), mcpb->getParents()[0]->getPDG() );
							ret = FindParton(mcpb);
						}
					}
					// for checking some strangge cases
					if (ret == 0 && abs(mcp->getPDG()) <= 6)
					{
						for (unsigned int i = 0; i < mcpa->getParents().size(); i++)
						{
							MCParticle *mcpb = mcpa->getParents()[i];
							printf("%5d vs. %5d, ", mcp->getPDG(), mcpb->getPDG());
						}
						printf("\n");
					}
				}
				else
				{
					// itself is a hadron
					// first; get all partons from the grandparents
					map<MCParticle *, TVector3> partonmap;
					partonmap.clear();
					map<MCParticle *, TVector3>::iterator it;
					for (unsigned int i = 0; i < mcpa->getParents().size(); i++)
					{
						MCParticle *mcpb = mcpa->getParents()[i];
						MCParticle *mcpo = NULL; //FindParton(mcpb);
						TVector3 v(mcpb->getMomentum());
						if (mcpo != NULL)
						{
							it = partonmap.find(mcpo);
							if (it != partonmap.end())
							{
								v += it->second; // why add?
												 //printf("the mother of %6d is %6d\n", (int)mcpb->getPDG(), (int)mcpo->getPDG());
							}
							partonmap[mcpo] = v;
						}
					}
					// then take the nearest one as the favorite parton
					double nearestangle = 2 * 3.14159;
					int selectedparton = -1, pdgid = 0;
					int pSelected = 0;
					int j;
					TVector3 v2(mcp->getMomentum());
					for (j = 0, it = partonmap.begin(); it != partonmap.end(); it++, j++)
					{
						if (it->first)
						{
							double angle = v2.Angle(it->second);
							if (nearestangle > angle)
							{
								nearestangle = angle;
								selectedparton = j;
								pSelected = (it->first)->getPDG();
								pdgid = (it->first)->getPDG();
							}
						}
					}
					ret = pSelected;
					if (ret == 0)
					{
						//cout << partonmap.size() << " partons obtained." << endl;
						//cout << "Parton " << selectedparton << " selected. "<< nearestangle  << endl;
					}
					FreeAll(partonmap);
				}
			}
			else
			{
				ret = FindParton(mcpa);
				if (0 && ret == 0)
				{
					for (unsigned int k = 0; k < mcpa->getParents().size(); k++)
					{
						MCParticle *part = mcpa->getParents()[k];
						printf("parent of %10d is %10d (%3d)\n", mcpa->getPDG(), part->getPDG(), (int)part->getParents().size());
						for (unsigned int l = 0; l < part->getParents().size(); l++)
						{
							printf("parent %2d is %10d\n", l, part->getParents()[l]->getPDG());
						}
					}
				}
			}
		}
	}
	else // More than one parent
	{
		if (abs(mcp->getPDG()) <= 6)
		{
			int mcpo2 = 0;
			for (unsigned int i = 0; i < mcp->getParents().size(); i++)
			{
				MCParticle *mcpa = mcp->getParents()[i];
				if (mcp->getParents().size() == 2 && abs(mcpa->getPDG()) == 11)
				{
					mcpo2 = mcp->getPDG();
					break;
				}
				else if (mcp->getPDG() == mcp->getPDG())
				{
					int mcpo = FindParton(mcpa);
					mcpo2 = mcpo;
					break;
				}
			}
			ret = mcpo2;
		}
		else
		{
			int mcpo2 = 0;
			for (unsigned int i = 0; i < mcp->getParents().size(); i++)
			{
				MCParticle *mcpa = mcp->getParents()[i];
				int mcpo = FindParton(mcpa);
				if (mcpo2 && mcpo != mcpo2)
				{
					cout << mcp->getPDG() << " has Multiple parents other than intermidiates with different parton origin!" << endl;
				}
				mcpo2 = mcpo;
			}
			//
			ret = mcpo2;
		}
	}
	//
	return ret;
}
