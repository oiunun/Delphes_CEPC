#include "NTupleHelper.h"
#include "Utilities.h"
const int NTupleHelper::maxsize =  100;
// ********************************************************************
//    EVENT INFORMATION
// ********************************************************************
void
NTupleHelper::fillEvent( LCEvent * evtP){
	fillDouble("Run",        evtP->getRunNumber()   );
	fillDouble("Event",      evtP->getEventNumber() );
}
//
void
NTupleHelper::fillExtras( vector<FSParticle*> PFOs, string cat, const int type){

	const double eps =1e-6;
	sort(PFOs.begin(), PFOs.end(), Sort_PFOs_E);
	vector<TLorentzVector> p4list;
	vector<double>         charge, vecVD, vecVZ;	
	//
	for (unsigned int j = 0; j < PFOs.size(); j++){
		FSParticle* fsp = PFOs[j];
		//
		double nsigd    = -99999;
		double nsigz    = -99999;
		if(fabs(fsp->charge())>0.1){
			TrackVec   tckvec = (fsp->pfo())->getTracks();
			if ( tckvec.size()>0){
				double d0       = tckvec[0]->getD0();
				double z0       = tckvec[0]->getZ0();
				double deltad0  = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
				double deltaz0  = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
				nsigd    = d0/(deltad0+eps);
				nsigz    = z0/(deltaz0+eps);
			}
		}
		//
		p4list.push_back(fsp->rawFourMomentum());
		charge.push_back(fsp->charge());
		vecVD.push_back(nsigd);
		vecVZ.push_back(nsigz);
	}
	//
	string index = "elec_idx";
	if      ( type == 11 ) index = "elec_idx";
	else if ( type == 13 ) index = "muon_idx";
	else if ( type == 22 ) index ="gamma_idx";
	else {
		cout<<"only for photon, muon or electron, but it is  "<< type<<endl;
		exit(1);
	}
	//
	fill4Momentum( index, cat , p4list,   PFOs.size());
	//
	if( type != 22 ){
		fillArray  ( concatName(cat, "charge"),  index, charge, PFOs.size());
		fillArray  ( concatName(cat,  "d0sig"),  index,  vecVD, PFOs.size());
		fillArray  ( concatName(cat,  "z0sig"),  index,  vecVZ, PFOs.size());
	}
	//
	FreeAll(p4list);
	FreeAll(charge);
	FreeAll( vecVD);
	FreeAll( vecVZ);
}
//
double NTupleHelper::JetCharge(ReconstructedParticle* jet, const double kappa){
		ReconstructedParticleVec vPart = jet->getParticles();
		int nPart = vPart.size(), id = -1 ; 
		TVector3 MOMENTUM  = TVector3(jet->getMomentum());
		double P           = 1; 
		if( id == 0 ) P    = jet->getEnergy();
		if( id == 1 ) P    = MOMENTUM.Pt();
		if( id == 2 ) P    = MOMENTUM.Mag();
		if(nPart>0){
			double charge = 0 , wate = 0 ;
			for (int k=0 ; k < nPart ; k++ ) {
				TVector3 momentum   = TVector3(vPart[k]->getMomentum());
				double c            = vPart[k]->getCharge();
				//cout<<c<<endl;
				if( fabs(c)<0.1 || fabs(c)>1.1)     continue;	
				double   p          = 0.0; 
				if( id <  0 ) p     = TMath::Cos(momentum.Angle(MOMENTUM));
				if( id == 0 ) p     = vPart[k]->getEnergy();
				if( id == 1 ) p     = momentum.Pt()/P;
				if( id == 2 ) p     = momentum.Mag();
				//if( fabs(p)<1.0 ) continue;	
				//
				double w  = pow( fabs(p), kappa );
				charge   += ( c*w );
				wate     += (  w  );
			}

			//cout << "wate/kappa = " << wate << "/" << kappa << endl;
			if( id <  0 ) charge = charge;//wate; 
			if( id == 0 ) charge = charge;//wate; 
			if( id == 1 ) charge = charge;//wate; 
			if( id == 2 ) charge = charge;//wate; 
			return charge;
		}
		return 9999.0;
}
//
void 
NTupleHelper::fillJet(FSParticle * fsp, int index, string tag , int m_full, int m_ft, double kappa){
	ReconstructedParticle * recPart = fsp->pfo();
	ReconstructedParticleVec vPart  = recPart->getParticles();
	//
	int ntrk = 0;
	int nclu = 0;
	int nPFO = vPart.size();
	double d0=0., z0=0., deltad0=0., deltaz0=0., nsigd0=0., nsigz0=0.;
	if( nPFO>0 && (m_full>0 || m_ft >0)){
		for (int k=0 ; k < nPFO ; k++ ) {
			TrackVec                 tckvec = vPart[k]->getTracks();
			ClusterVec               cluvec = vPart[k]->getClusters();
			ntrk   += tckvec.size();
			nclu   += cluvec.size();
			if (tckvec.size() > 0) {
				for (unsigned int j=0 ; j < tckvec.size() ; j++ ) {
					d0      += tckvec[j]->getD0();
					z0      += tckvec[j]->getZ0();
					deltad0  = TMath::Sqrt(tckvec[j]->getCovMatrix()[0]);
					deltaz0  = TMath::Sqrt(tckvec[j]->getCovMatrix()[9]);
					nsigd0  += tckvec[j]->getD0()/deltad0;
					nsigz0  += tckvec[j]->getZ0()/deltaz0;
				}
			}
		}
		d0      /= ntrk;
		z0      /= ntrk;
		nsigd0  /= ntrk;
		nsigz0  /= ntrk;
	}
	//
	if(0) printf("Ntrk = %3d, Nclu = %3d, Npar = %3d\n", ntrk, nclu, nPFO);
	//
	double    energy    = recPart->getEnergy();
	double  charge      = recPart->getCharge();
	//double    charge     = JetCharge(recPart, kappa);
	double   mass       = recPart->getMass();
	TVector3 momentum   = TVector3(recPart->getMomentum());
	HepLorentzVector p4mom =HepLorentzVector (momentum.X(), momentum.Y(), momentum.Z(),energy);
	double   momentumM  = momentum.Mag();
	//double   rapidity   = 0.5*TMath::Log((energy+momentum.Pz())/(energy-momentum.Pz()));
	double   rapidity   = momentum.PseudoRapidity();
	double   cosTheta   = momentum.CosTheta();
	//double   sphericity = fsp->sphericity(); 
	//
	if( m_full > 0 || m_ft ){
		fillDouble(concatName(tag,"VtxR"     ,  index), d0              );
		fillDouble(concatName(tag,"VtxZ"     ,  index), z0              );
		fillDouble(concatName(tag,"VtxSigR"  ,  index), nsigd0          );
		fillDouble(concatName(tag,"VtxSigZ"  ,  index), nsigz0          );
		fillDouble(concatName(tag,"Btag"     ,  index), fsp->btag()     );
		fillDouble(concatName(tag,"Ctag"     ,  index), fsp->ctag()     );
		fillDouble(concatName(tag,"Cat"      ,  index), fsp->cat()      );
	}
	fillDouble(concatName(tag,"ntrk"      ,  index), ntrk            );
	fillDouble(concatName(tag,"nclu"      ,  index), nclu            );
	fillDouble(concatName(tag,"charge"    ,  index), charge          );
	fillDouble(concatName(tag,"nPFO"      ,  index), nPFO            );
	fillDouble(concatName(tag,"mass"      ,  index), mass            );
	fillDouble(concatName(tag,"En"        ,  index), energy          );
	fillDouble(concatName(tag,"Px"        ,  index), momentum.X()    );
	fillDouble(concatName(tag,"Py"        ,  index), momentum.Y()    );
	fillDouble(concatName(tag,"Pz"        ,  index), momentum.Z()    );
	fillDouble(concatName(tag,"Pt"        ,  index), momentum.Pt()   );
	fillDouble(concatName(tag,"Ptot"      ,  index), momentumM       );
	fillDouble(concatName(tag,"Rapidity"  ,  index), rapidity        );
	fillDouble(concatName(tag,"cosTheta"  ,  index), cosTheta        );
	//fillDouble(concatName(tag,"Sphericity",  index), sphericity      );
	double pdgid=0;
	HepLorentzVector p4(0,0,0,0);
	MCParticle * mcp= fsp->mcp();
	if( mcp ) {
		pdgid= mcp->getPDG();
		p4=HepLorentzVector (mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], mcp->getEnergy() );
	} 
	fillDouble(concatName(tag,"PDGID"      ,  index),   pdgid        );
	fillDouble(concatName(tag,"McPx"       ,  index),   p4.px()      );
	fillDouble(concatName(tag,"McPy"       ,  index),   p4.py()      );
	fillDouble(concatName(tag,"McPz"       ,  index),   p4.pz()      );
	fillDouble(concatName(tag,"McEn"       ,  index),   p4.e ()      );
	fillDouble(concatName(tag,"AngleRecMc" ,  index),   p4.angle(p4mom.v()));
}
void 
NTupleHelper::fillTracks( vector<Track*> tckvec ){

	vector<double> VecD0, VecZ0, VecD0Sig, VecZ0Sig, VecCos,VecP;
	double  d0=9990.,z0=9990.,deltad0=0.,deltaz0=0.,nsigd0=9990.,nsigz0=9990.,cos=999,pp=-999;
	for (unsigned int j=0 ; j < tckvec.size() ; j++ ) {
		d0       = tckvec[j]->getD0();
		z0       = tckvec[j]->getZ0();
		deltad0  = TMath::Sqrt(tckvec[j]->getCovMatrix()[0]);
		deltaz0  = TMath::Sqrt(tckvec[j]->getCovMatrix()[9]);
		nsigd0   = tckvec[j]->getD0()/deltad0;
		nsigz0   = tckvec[j]->getZ0()/deltaz0;
		cos      = tckvec[j]->getTanLambda();
		cos      = 1/pow(1+pow(cos,2),0.5)*(cos/fabs(cos));
		VecD0 .push_back(d0);
		VecZ0 .push_back(z0);
		VecCos.push_back(cos);
		//VecP  .push_back(pp);
	}
	
	fillArray ("trkD0"    , "idx_trk", VecD0    , tckvec.size());
	fillArray ("trkZ0"    , "idx_trk", VecZ0    , tckvec.size());
	fillArray ("trkCos"   , "idx_trk", VecCos   , tckvec.size());
}

void 
NTupleHelper::fillPFO(ReconstructedParticle * recPart, FSParticle * fsp, int index, string tag, int m_full ){
	TrackVec   tckvec = recPart->getTracks();
	ClusterVec cluvec = recPart->getClusters();
	int ntrk = tckvec.size();
	int nclu = cluvec.size();
	double  d0=9990.,z0=9990.,deltad0=0.,deltaz0=0.,nsigd0=9990.,nsigz0=9990.;
	if (m_full > 0 && ntrk > 0) {
		d0=0., z0=0., nsigd0=0., nsigz0=0.;
		for (unsigned int j=0 ; j < tckvec.size() ; j++ ) {
			d0      += tckvec[j]->getD0();
			z0      += tckvec[j]->getZ0();
			deltad0  = TMath::Sqrt(tckvec[j]->getCovMatrix()[0]);
			deltaz0  = TMath::Sqrt(tckvec[j]->getCovMatrix()[9]);
			nsigd0  += tckvec[j]->getD0()/deltad0;
			nsigz0  += tckvec[j]->getZ0()/deltaz0;
		}
		d0      /= ntrk;
		z0      /= ntrk;
		nsigd0  /= ntrk;
		nsigz0  /= ntrk;
	}
	//
	double   leptonType = fsp->leptonType(); 
	double   energy     = recPart->getEnergy();
	double   charge     = recPart->getCharge();
	double   mass       = recPart->getMass();
	TVector3 momentum   = TVector3(recPart->getMomentum());
	double   momentumM  = momentum.Mag();
	HepLorentzVector p4(momentum.X(), momentum.Y(), momentum.Z(),energy);
	//double   rapidity   = 0.5*TMath::Log((energy+momentum.Pz())/(energy-momentum.Pz()));
	double   rapidity   = momentum.PseudoRapidity();
	double   cosTheta   = momentum.CosTheta();
	double   ecalEnergy = 0;
	double   hcalEnergy = 0;
	double   yokeEnergy = 0;
	double   totalCalEnergy = 0;
	int nHits = 0;
	if ( m_full>0  && energy > 0) {
		std::vector<lcio::Cluster*> clusters = recPart->getClusters();
		for (std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster) {
			ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
			hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
			yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
			ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
			hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
			CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
			nHits += calHits.size();
		}
		totalCalEnergy = ecalEnergy + hcalEnergy;
	}
	//
	if ( m_full>0 ) {
		fillDouble(concatName(tag,"ntrk"     ,  index), ntrk);
		fillDouble(concatName(tag,"nclu"     ,  index), nclu);
		fillDouble(concatName(tag,"VtxR"     ,  index), d0);
		fillDouble(concatName(tag,"VtxZ"     ,  index), z0);
		fillDouble(concatName(tag,"ecal"     ,  index), ecalEnergy  );
		fillDouble(concatName(tag,"hcal"     ,  index), hcalEnergy  );
		fillDouble(concatName(tag,"nHits"    ,  index), nHits);
		fillDouble(concatName(tag,"totCalEn" ,  index), totalCalEnergy);
	}
	fillDouble(concatName(tag,"charge"   ,  index), charge);
	fillDouble(concatName(tag,"mass"     ,  index), mass  );
	fillDouble(concatName(tag,"En"       ,  index), energy);
	fillDouble(concatName(tag,"Px"       ,  index), momentum.X());
	fillDouble(concatName(tag,"Py"       ,  index), momentum.Y());
	fillDouble(concatName(tag,"Pz"       ,  index), momentum.Z());
	fillDouble(concatName(tag,"Pt"       ,  index), momentum.Pt());
	fillDouble(concatName(tag,"Ptot"     ,  index), momentumM   );
	fillDouble(concatName(tag,"LepType"  ,  index), leptonType);
	fillDouble(concatName(tag,"Rapidity" ,  index), rapidity    );
	fillDouble(concatName(tag,"cosTheta" ,  index), cosTheta    );
	double pdgid=0, mcten=0, mctpx=0, mctpy=0, mctpz=0;
	MCParticle * mcp= fsp->mcp();
	if( mcp ) {
		//cout<<"I am here to check mcp "<<endl;
		mcten  = mcp->getEnergy();
		TVector3  p3  = TVector3(mcp->getMomentum());
		mcten  = mcp->getEnergy();
		mctpx  = p3.X(); 
		mctpy  = p3.Y(); 
		mctpz  = p3.Z(); 
		HepLorentzVector p4mc(mctpx,mctpy,mctpz,mcten);
		pdgid= mcp->getPDG();
		fillDouble(concatName(tag,"PDGID"      ,  index),   pdgid  );
		fillDouble(concatName(tag,"MCTEN"      ,  index),   mcten  );
		fillDouble(concatName(tag,"MCTPX"      ,  index),   mctpx  );
		fillDouble(concatName(tag,"MCTPY"      ,  index),   mctpy  );
		fillDouble(concatName(tag,"MCTPZ"      ,  index),   mctpz  );
		fillDouble(concatName(tag,"AngleRecMc" ,  index),   p4.angle(p4mc.v()));
	}	
}
void 
NTupleHelper::fillPFOs(vector<ReconstructedParticle*> pfos, LCCollection *colMCTL){

	vector<double> VecEn, VecPp, VecPx, VecPy, VecPz, VecCs, VecCh, VecPdgid;
	for (unsigned int j = 0; j < pfos.size(); j++){
		ReconstructedParticle * pfo= pfos[j];
		double   energy     = pfo->getEnergy();
		double   charge     = pfo->getCharge();
		TVector3 momentum   = TVector3(pfo->getMomentum());
		double   momentumM  = momentum.Mag();
		double   cosTheta   = momentum.CosTheta();
		//
		VecEn.push_back(energy);
		VecPp.push_back(momentumM);
		VecPx.push_back(momentum.X());
		VecPy.push_back(momentum.Y());
		VecPz.push_back(momentum.Z());
		VecCs.push_back(cosTheta);
		VecCh.push_back(charge);
		double pdgid=0;
		if ( colMCTL ){
			LCRelationNavigator *navMCTL   = new LCRelationNavigator(colMCTL);
			LCObjectVec vecMCTL            = navMCTL->getRelatedToObjects(pfo);
			if(0) {
				printf("vecMCTL size = %3d\n", (int)vecMCTL.size());
			}
			if (vecMCTL.size() > 0) {
				MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL[0]);
				if(mcp) pdgid=mcp->getPDG();
			}
		delete navMCTL;
		}
		VecPdgid.push_back(pdgid);
	}
	//
	fillArray ("Pfo_Pdgid" , "idx_pfo", VecPdgid , pfos.size());
	fillArray ("Pfo_En"    , "idx_pfo", VecEn    , pfos.size());
	fillArray ("Pfo_Pp"    , "idx_pfo", VecPp    , pfos.size());
	fillArray ("Pfo_Px"    , "idx_pfo", VecPx    , pfos.size());
	fillArray ("Pfo_Py"    , "idx_pfo", VecPy    , pfos.size());
	fillArray ("Pfo_Pz"    , "idx_pfo", VecPz    , pfos.size());
	fillArray ("Pfo_Cos"   , "idx_pfo", VecCs    , pfos.size());
	fillArray ("Pfo_Charge", "idx_pfo", VecCh    , pfos.size());

}
void 
NTupleHelper::fillPFO(FSParticle * fsPart, int index, string tag){
	//
	double   energy     = fsPart->energy();
	double   charge     = fsPart->charge();
	TVector3 momentum   = TVector3((fsPart->rawFourMomentum()).Vect());
	double   momentumM  = momentum.Mag();
	double   rapidity   = 0.5*TMath::Log((energy+momentum.Pz())/(energy-momentum.Pz()));
	double   cosTheta   = momentum.CosTheta();
	//
	fillDouble(concatName(tag,"charge"   ,  index), charge);
	fillDouble(concatName(tag,"En"       ,  index), energy);
	fillDouble(concatName(tag,"Px"       ,  index), momentum.X());
	fillDouble(concatName(tag,"Py"       ,  index), momentum.Y());
	fillDouble(concatName(tag,"Pz"       ,  index), momentum.Z());
	fillDouble(concatName(tag,"Pt"       ,  index), momentum.Pt());
	fillDouble(concatName(tag,"Ptot"     ,  index), momentumM   );
	fillDouble(concatName(tag,"Rapidity" ,  index), rapidity    );
	fillDouble(concatName(tag,"cosTheta" ,  index), cosTheta    );
}

//
void
NTupleHelper::fill4Momentum(int index, string tag, const TLorentzVector& p){
	fillDouble(concatName(tag,"Px",index), p.Px());
	fillDouble(concatName(tag,"Py",index), p.Py());
	fillDouble(concatName(tag,"Pz",index), p.Pz());
	fillDouble(concatName(tag,"En",index), p.E() );
}
//
void
NTupleHelper::fill4Momentum(string tag, const TLorentzVector& p){
	fillDouble(concatName(tag,"Px",0), p.Px());
	fillDouble(concatName(tag,"Py",0), p.Py());
	fillDouble(concatName(tag,"Pz",0), p.Pz());
	fillDouble(concatName(tag,"En",0), p.E() );
}
//
void
NTupleHelper::fill4Momentum(int index, string tag){
	fillDouble(concatName(tag,"Px",index), 0.0);
	fillDouble(concatName(tag,"Py",index), 0.0);
	fillDouble(concatName(tag,"Pz",index), 0.0);
	fillDouble(concatName(tag,"En",index), 0.0);
}
//
// ********************************************************************
//    FOUR MOMENTA
// ********************************************************************
void
NTupleHelper::fill4Momentum(string index_name, string tag, const vector<TLorentzVector>& p, const int size){

	int size_loc=size;
	if(size > maxsize){
		//printf("too many entries %3d (maximum is 100), please check\n", size);
		//exit(1);
		size_loc=maxsize;
	}	
	if (m_bookingStage && !containsEntry(index_name)){
		char Name[20], Description[20];
		sprintf(Name,"%s",index_name.c_str());
		sprintf(Description,"%s/I",index_name.c_str());
		m_Tree->Branch(Name , &m_ntupleIntMap[index_name], Description);	
	}
	if (!m_bookingStage && !containsEntry(index_name)){
		cout << "NTUPLEHELPER:  Variable " << index_name << " has not been booked." << endl;
		exit(0);
	}
	m_IntMap   [index_name] = size_loc;
	//
	//
	string combname[4]={"Px","Py","Pz","En"}, name;
	for(int i=0; i<4; i++){
		name=concatName(tag,combname[i]);
		if (m_bookingStage && !containsEntry(name)){
			char Name[20], Description[20];
			sprintf(Name       ,"%s[%d]"  , name.c_str(), maxsize);
			sprintf(Description,"%s[%d]/D", name.c_str(), maxsize);
			m_ntupleArrayMap[name]= new double[maxsize];
			m_Tree->Branch(Name, m_ntupleArrayMap[name], Description);	
		}
		if (!m_bookingStage && !containsEntry(name)){
			cout << "NTUPLEHELPER:  Variable " << name << " has not been booked." << endl;
			exit(0);
		}
		m_arrayMap [name] = new double[maxsize];
		memset(m_arrayMap[name], 0, sizeof(double)*maxsize);
		m_arraySize[name] = size_loc;
		for(int j=0; j<size_loc; j++){
				if(i==0) m_arrayMap[name][j] = p[j].Px();
				if(i==1) m_arrayMap[name][j] = p[j].Py();
				if(i==2) m_arrayMap[name][j] = p[j].Pz();
				if(i==3) m_arrayMap[name][j] = p[j].E ();
		}
	}
}

// ********************************************************************
//    HELPER FUNCTIONS
// ********************************************************************
NTupleHelper::NTupleHelper(TTree* Tree, TLorentzVector p4){
	m_Tree         = Tree;
	p4psi          =   p4;
	m_bookingStage = true;
	FreeAll(m_doubleMap);
	FreeAll(m_IntMap   );
	FreeAll(m_arrayMap );
	FreeAll(m_arraySize);
	if( m_Tree == NULL )
		cout    << "ERROR:  null tree pointer -- "
			<< "check for duplicate final states in configuration"
			<< endl; assert( m_Tree != NULL );
}

NTupleHelper::~NTupleHelper(){
	m_Tree->Write();
	for (map<string,double*>::iterator mapItr = m_ntupleArrayMap.begin();
			mapItr != m_ntupleArrayMap.end(); mapItr++){
		delete []mapItr->second;
	}
	delete m_Tree;
}

void
NTupleHelper::write(){
	
	for (map<string,double>::iterator mapItr = m_doubleMap.begin();
			mapItr != m_doubleMap.end(); mapItr++){
		m_ntupleDoubleMap[mapItr->first] = mapItr->second;
	}
	
	for (map<string,int>::iterator mapItr = m_IntMap.begin();
			mapItr != m_IntMap.end(); mapItr++){
		m_ntupleIntMap[mapItr->first] = mapItr->second;
	}
	
	for (map<string,double*>::iterator mapItr = m_arrayMap.begin();
			mapItr != m_arrayMap.end(); mapItr++){
		for( int i=0; i< maxsize; i++)
		{
			m_ntupleArrayMap[mapItr->first][i]=mapItr->second[i];
		}
	}

	m_Tree->Fill();

	for (map<string,double*>::iterator mapItr = m_arrayMap.begin();
			mapItr != m_arrayMap.end(); mapItr++){
		delete[] mapItr->second;
	}

	m_bookingStage = false;
	
}


void NTupleHelper::fillDouble(string name, double value){
	if (m_bookingStage && !containsEntry(name)){
		char Name[20], Description[20];
		sprintf(Name,"%s",name.c_str());
		sprintf(Description,"%s/D",name.c_str());
		m_Tree->Branch(Name , &m_ntupleDoubleMap[name], Description);	
	}
	if (!m_bookingStage && !containsEntry(name)){
		cout << "NTUPLEHELPER:  Variable " << name << " has not been booked." << endl;
		exit(0);
	}
	m_doubleMap[name] = value;
}

void NTupleHelper::fillLong(string name, int value){
	if (m_bookingStage && !containsEntry(name)){
		char Name[20], Description[20];
		sprintf(Name,"%s",name.c_str());
		sprintf(Description,"%s/I",name.c_str());
		m_Tree->Branch(Name   , &m_ntupleIntMap[name], Description);	
	}
	if (!m_bookingStage && !containsEntry(name)){
		cout << "NTUPLEHELPER:  Variable " << name << " has not been booked." << endl;
		exit(0);
	}
	m_IntMap[name] = value;
}


void 
NTupleHelper::fillArray(string name, string index_name, double* value, int size){
	
	int size_loc=size;
	if(size >maxsize){
		//printf("too many entries %3d (maximum is 30), please check\n", size);
		//exit(1);
		size_loc=maxsize;
	}	
	if (m_bookingStage && !containsEntry(index_name)){
		char Name[20], Description[20];
		sprintf(Name,"%s",index_name.c_str());
		sprintf(Description,"%s/I",index_name.c_str());
		m_Tree->Branch(Name, &m_ntupleIntMap[index_name], Description);
	}
	
	m_IntMap[index_name]=size_loc;
	m_arraySize[name   ]=size_loc;

	if (m_bookingStage && !containsEntry(name)){
		char Name[20], Description[20];
		sprintf(Name       ,"%s[%d]"  , name.c_str(), maxsize);
		sprintf(Description,"%s[%d]/D", name.c_str(), maxsize);
		m_ntupleArrayMap[name]= new double[maxsize];
		m_Tree->Branch(Name, m_ntupleArrayMap[name], Description);	
	}
	
	if (!m_bookingStage && !containsEntry(name)){
		cout << "NTUPLEHELPER:  Variable " << name << " has not been booked." << endl;
		exit(0);
	}
	
	if (!m_bookingStage && !containsEntry(index_name)){
		cout << "NTUPLEHELPER:  Variable " << index_name << " has not been booked." << endl;
		exit(0);
	}
        
	m_arrayMap[name]= new double[maxsize];
	memset(m_arrayMap[name], 0, sizeof(double)*maxsize);
	for(int i=0; i<size_loc; i++){
		m_arrayMap[name][i] = *(value+i);
	}
}

void 
NTupleHelper::fillArray(string name, string index_name, vector<double> value, int size){

	int size_loc=size;
	if(size >maxsize){
		//printf("too many entries %3d (maximum is 30), please check\n", size);
		//exit(1);
		size_loc=maxsize;
	}	
	if (m_bookingStage && !containsEntry(index_name)){
		char Name[20], Description[20];
		sprintf(Name,"%s",index_name.c_str());
		sprintf(Description,"%s/I",index_name.c_str());
		m_Tree->Branch(Name, &m_ntupleIntMap[index_name], Description);	
	}
	
	m_IntMap[index_name]=size_loc;
	m_arraySize[name   ]=size_loc;

	if (m_bookingStage && !containsEntry(name)){
		char Name[20], Description[20];
		sprintf(Name       ,"%s[%d]"  , name.c_str(), maxsize);
		sprintf(Description,"%s[%d]/D", name.c_str(), maxsize);
		m_ntupleArrayMap[name]= new double[maxsize];
		m_Tree->Branch(Name, m_ntupleArrayMap[name], Description);	
	}
	
	if (!m_bookingStage && !containsEntry(name)){
		cout << "NTUPLEHELPER:  Variable " << name << " has not been booked." << endl;
		exit(0);
	}
	
	if (!m_bookingStage && !containsEntry(index_name)){
		cout << "NTUPLEHELPER:  Variable " << index_name << " has not been booked." << endl;
		exit(0);
	}

	m_arrayMap[name]= new double[maxsize];
	memset(m_arrayMap[name], 0, sizeof(double)*maxsize);
	for(int i=0; i<size_loc; i++){
		m_arrayMap[name][i] = value[i];
	}
}

bool
NTupleHelper::containsEntry(string name){
	map<string,double>::iterator mapItr1 = m_doubleMap.find(name);
	if (mapItr1 != m_doubleMap.end()) return true;
	
	map<string,int>::iterator mapItr2 = m_IntMap.find(name);
	if (mapItr2 != m_IntMap.end()   ) return true;
	
	map<string,double*>::iterator mapItr3 = m_arrayMap.find(name);
	if (mapItr3 != m_arrayMap.end() ) return true;
	
	return false;
}

string 
NTupleHelper::concatName(string tag, string base, int index){
	stringstream name;
	name << tag;
	name << base;
	if( index>0){
		name << "P";
		name << index;
	}
	return name.str();
}

string 
NTupleHelper::concatName(string tag, string base){
	stringstream name;
	name << tag;
	name << base;
	return name.str();
}

void //
NTupleHelper::fillMCTruth(MCTruthHelper* mc ){
	const int limit = 999;
	vector<MCParticleImpl*> mcPList    = mc->finalParticleList();
	//
	int    m_numParticle=0;
	double mc_pdgid[1000], mc_motheridx[1000], mc_trkidx[1000];
	double mc_px[1000], mc_py[1000], mc_pz[1000], mc_pe[1000];
	for (unsigned int i=0; i < mcPList.size(); i++)
	{
		int nparents   = mcPList[i]->getParents().size();
		int pdgid      = mcPList[i]->getPDG();
		int id         = mcPList[i]->id();
		Int_t mother = 0, mid=0;
		if (nparents > 0) {
			MCParticle *Parent = mcPList[i]->getParents()[0];
			mother             = Parent->getPDG();
			mid                = Parent->id();
		}
		if (mc->hasMothers( mcPList[i],111)) continue;
		if (mc->hasMothers( mcPList[i],221)) continue;
		if (abs(pdgid)==12||abs(pdgid)==14||abs(pdgid)==16) continue;
		   
		TLorentzVector p4p   = TLorentzVector(mcPList[i]->getMomentum(),mcPList[i]->getEnergy());
		mc_trkidx   [m_numParticle] = id ;
		mc_motheridx[m_numParticle] = mid ;
		mc_pdgid    [m_numParticle] = pdgid  ;
		mc_px       [m_numParticle] = p4p.X();
		mc_py       [m_numParticle] = p4p.Y();
		mc_pz       [m_numParticle] = p4p.Z();
		mc_pe       [m_numParticle] = p4p.T();
		m_numParticle++;
		if(m_numParticle>limit )
		{
			cout<<"m_numParticle >  : "<< limit << " " <<m_numParticle<<endl;
			break;
		}
	}
	//
	//fillArray ("trkidx"    , "idx_mc", (double*)mc_trkidx    , m_numParticle);
	//fillArray ("motheridx" , "idx_mc", (double*)mc_motheridx , m_numParticle);
	fillArray ("mc_pdgid"  , "idx_mc", (double*)mc_pdgid     , m_numParticle);
	fillArray ("mc_En"     , "idx_mc", (double*)mc_pe        , m_numParticle);
	fillArray ("mc_Px"     , "idx_mc", (double*)mc_px        , m_numParticle);
	fillArray ("mc_Py"     , "idx_mc", (double*)mc_py        , m_numParticle);
	fillArray ("mc_Pz"     , "idx_mc", (double*)mc_pz        , m_numParticle);
}

void // this function is not ready due to the complication of whizard MC topology  
NTupleHelper::fillMCTruthDT(MCTruthHelper* mc, FSInfo* fs ){
	//printf("\n\n"); 
	if (fs){
		const int limit = 999;
		vector<MCParticleImpl*> mcPList    = mc->AllParticleList();
		//
		int    m_numParticle=0;
		double mc_pdgid[1000], mc_motheridx[1000], mc_trkidx[1000];
		double mc_px[1000], mc_py[1000], mc_pz[1000], mc_pe[1000];
		for (unsigned int i=0; i < mcPList.size(); i++)
		{
			//int isCreatedInSimulation         = mcPList[i]->isCreatedInSimulation()       ;
			//int isBackscatter                 = mcPList[i]->isBackscatter()               ; 
			//int isOverlay                     = mcPList[i]->isOverlay()                   ;
			//int vertexIsNotEndpointsOfParent  = mcPList[i]->vertexIsNotEndpointOfParent() ; 
			int nparents   = mcPList[i]->getParents().size();
			int pdgid      = mcPList[i]->getPDG();
			int id         = mcPList[i]->id();
			Int_t mother = 0, mid=0;
			if (nparents > 0) {
				MCParticle *Parent = mcPList[i]->getParents()[0];
				mother             = Parent->getPDG();
				mid                = Parent->id();
			}
			//if (    mcPList[i]->isCreatedInSimulation()        ) continue;
			//if (    mcPList[i]->isBackscatter()                ) continue;
			//if (    mcPList[i]->isOverlay()                    ) continue;
			//if (    mcPList[i]->vertexIsNotEndpointOfParent()  ) continue;
			if (mc->hasMothers( mcPList[i],111) ) continue;
			if (mc->hasMothers( mcPList[i],221) ) continue;
			/*
			printf("No. %3d, pdgid = %8d %8d, status = %2d %2d %2d %2d, mass=%9.4f\n", i, pdgid,mother,
					isCreatedInSimulation,      
					isBackscatter         ,     
					isOverlay              ,    
					vertexIsNotEndpointsOfParent,
					mcPList[i]->getMass());
			*/
			TLorentzVector p4p   = TLorentzVector(mcPList[i]->getMomentum(),mcPList[i]->getEnergy());
			mc_trkidx   [m_numParticle] = id ;
			mc_motheridx[m_numParticle] = mid ;
			mc_pdgid    [m_numParticle] = pdgid  ;
			mc_px       [m_numParticle] = p4p.X();
			mc_py       [m_numParticle] = p4p.Y();
			mc_pz       [m_numParticle] = p4p.Z();
			mc_pe       [m_numParticle] = p4p.T();
			m_numParticle++;
			if(m_numParticle>limit )
			{
				cout<<"m_numParticle >  : "<< limit << " " <<m_numParticle<<endl;
				break;
			}
		}
		//
		//fillArray ("trkidx"    , "indexmc", (double*)mc_trkidx    , m_numParticle);
		//fillArray ("motheridx" , "indexmc", (double*)mc_motheridx , m_numParticle);
		fillArray ("mc_pdgid"  , "indexmc", (double*)mc_pdgid     , m_numParticle);
		fillArray ("mc_Px"     , "indexmc", (double*)mc_px        , m_numParticle);
		fillArray ("mc_Py"     , "indexmc", (double*)mc_py        , m_numParticle);
		fillArray ("mc_Pz"     , "indexmc", (double*)mc_pz        , m_numParticle);
		fillArray ("mc_En"     , "indexmc", (double*)mc_pe        , m_numParticle);
	}
}
