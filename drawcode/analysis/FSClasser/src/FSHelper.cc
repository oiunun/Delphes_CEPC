#include "FSHelper.h"
//**********************************************
//
//   Missing particle
//
//**********************************************
FSParticle::FSParticle(string name, bool missed):
	m_name(name), m_missed(missed) 
{
	m_mass = Mass(m_name);
	m_rawFourMomentum = TLorentzVector(0,0,0,0);
}
//**********************************************
//**********************************************
FSParticle::FSParticle(
		ReconstructedParticle* pfo, LCCollection *colMCTL, 
		vector<MCParticleImpl*> partonList, string name,
		double btag, double ctag, double bctag, double cat, double flavor):
	m_name(name), m_missed(false), m_fast(false) 
{
	m_pfo                   = pfo;
	m_mass                  = Mass(m_name);
	m_rawFourMomentum       = fourMomentum(pfo);
	m_charge                = pfo->getCharge();
	m_type                  = pfo->getType();
	m_btag                  = btag;
	m_ctag                  = ctag;
	m_cat                   = cat ;
	m_bctag                 = bctag;
	m_flavor                = flavor;
	m_leptonType            = 0.0;
	m_Energy                = m_rawFourMomentum.E();
	m_pT                    = m_rawFourMomentum.Pt();
	m_pZ                    = m_rawFourMomentum.Pz();
	m_CosTheta              = m_rawFourMomentum.CosTheta();
	m_Rapidity              = m_rawFourMomentum.Rapidity();

	m_pdgid                 = PdgCode(m_name);
	m_mcp                   = NULL;
	if ( m_type!=4 ){ // not a jet
		TrackVec trkvec    = pfo->getTracks();
		ClusterVec cltvec  = pfo->getClusters();
		ReconstructedParticleVec parvec=  pfo->getParticles();
		Int_t ntrk    = trkvec.size();
		Int_t nclu    = cltvec.size();
		Int_t npar    = parvec.size();
		if(0) {
			printf("%10s: Ntrk = %3d, Nclu = %3d, Npart = %3d\n", m_name.data(), ntrk, nclu, npar);
			if (npar > 0) {
				for( int i=0;i<npar; i++){
					printf(" par # = %3d, %10d\n", i, parvec[i]->id() );
				}	  
				for( int i=0;i<ntrk; i++){
					printf(" trk # = %3d, %10d\n", i, trkvec[i]->id() );
				}	  
				for( int i=0;i<nclu; i++){
					printf(" clu # = %3d, %10d\n", i, cltvec[i]->id() );
				}	  
			}
		}
		if (npar > 0) {
			for( int i=0;i<npar; i++){
				m_particleId.push_back(parvec[i]->id());
			}	  
		}
		if (ntrk > 0) {
			for( int i=0;i<ntrk; i++){
				m_trackId.push_back(trkvec[i]->id());
			}	  
		}
		if (nclu > 0) {
			for( int i=0;i<nclu; i++){
				m_showerId.push_back(cltvec[i]->id());
			}	  
		}
		//
		if ( colMCTL ){
			LCRelationNavigator *navMCTL   = new LCRelationNavigator(colMCTL);
			LCObjectVec vecMCTL            = navMCTL->getRelatedToObjects(pfo);
			if(0) {
				printf("vecMCTL size = %3d\n", (int)vecMCTL.size());
			}
			//
			if (vecMCTL.size() > 0) {
				int pdgid=0;
				for(unsigned int k=0; k< vecMCTL.size(); k++){
					MCParticle* mcp = dynamic_cast<MCParticle *>(vecMCTL[k]);
					pdgid=mcp->getPDG();
					if (pdgid==m_pdgid){
						m_mcp = mcp;	
						break;
					}
				}
				if( m_pdgid != pdgid ) printf("recon object is %8d, mctruth is %8d\n", pdgid, m_pdgid);
			}
			//
			delete navMCTL;
			//if ( abs(m_pdgid)==11 || abs(m_pdgid)==13 ) m_leptonType=FindLeptonParent(m_mcp);
		}
	}else{
		m_pdgid = 0;
		m_mass  = pfo->getMass();
		ReconstructedParticleVec vPart = pfo->getParticles();
		int nPart= vPart.size();
		if(0) printf("Npar of jet = %3d\n", nPart);
		if(nPart>0){
			if (nPart > 0) {
				for( int i=0;i<nPart; i++){
					m_particleId.push_back(vPart[i]->id());
				}	  
			}
			for (int k=0 ; k < nPart ; k++ ) {
				TrackVec   tckvec               = vPart[k]->getTracks();
				ClusterVec cluvec               = vPart[k]->getClusters();
				ReconstructedParticleVec parvec = vPart[k]->getParticles();
				int npar   = parvec.size();
				int ntrk   = tckvec.size();
				int nclu   = cluvec.size();
				if (0) printf("ntrk = %3d, nclu = %3d, npar = %3d\n", ntrk, nclu, npar);
				if (ntrk > 0) {
					for(int it=0; it<ntrk; it++) 
						m_trackId.push_back(tckvec[it]->id());
				}
				if (nclu > 0) {
					for(int it=0; it<nclu; it++) 
						m_showerId.push_back(cluvec[it]->id());
				}
				if (npar > 0) {
					for( int i=0;i<npar; i++)
						m_particleId.push_back(parvec[i]->id());
				}
			}
			//
			if ( colMCTL ) {
				// matching to mc parton
				if ( 0 ) printf("Search for the parton monther of %2d reconstructed objects\n",nPart);
				map<MCParticle *, int > parentmap;
				for (int k=0 ; k < nPart ; k++ ) {
					ReconstructedParticle* part    = vPart[k];
					LCRelationNavigator *navMCTL   = new LCRelationNavigator(colMCTL);
					LCObjectVec vecMCTL            = navMCTL->getRelatedToObjects(part);
					if ( 0 ) printf("Size of parton list is %3d\n", (int)partonList.size());
					if ( 0 ) printf("matched MC objects of %2d is %2d\n",k, (int)vecMCTL.size() );
					if (vecMCTL.size() > 0) {
						for ( unsigned int j=0; j<vecMCTL.size(); j++ ) {
							MCParticle *mcp = dynamic_cast<MCParticle *>(vecMCTL[j]);
							//find parent parton
							MCParticle *mcparent = FindParton(mcp);
							if( mcparent != NULL ){
								map<MCParticle *, int>::iterator it;
								it = parentmap.find(mcparent);
								if( it != parentmap.end() ) 
									parentmap[mcparent] += 1;
								else
									parentmap[mcparent]  = 1;
							}
						}
					}
					delete navMCTL;
				}
				//
				if( 0 && parentmap.size()==0 ){
				       	printf("there is %3d particles in the jet, but failed to find the parton mother\n", nPart);
					for (int k=0 ; k < nPart ; k++ ) {
						ReconstructedParticle* part    = vPart[k];
						LCRelationNavigator *navMCTL   = new LCRelationNavigator(colMCTL);
						LCObjectVec vecMCTL            = navMCTL->getRelatedToObjects(part);
						if (vecMCTL.size() > 0) {
							for ( unsigned int j=0; j<vecMCTL.size(); j++ ) {
								MCParticle *mcp = dynamic_cast<MCParticle *>(vecMCTL[j]);
								printf("daughter %2d is %10d\n", k, mcp->getPDG() );
							}
						}
					}
				}
				map<MCParticle *, int>::iterator it;
				int ic=-99;
				for(it = parentmap.begin();it != parentmap.end();it++){
					if ( it->second > ic && it->first->getPDG()!=0  ){
						m_mcp   = it->first;
						m_pdgid = it->first->getPDG();
						ic      = it->second;	
					}
				}
				FreeAll(parentmap);

				if( 0 && m_mcp ) printf("matched quark is %3d\n", m_pdgid);
			}
			//
		}
	}
	//
	double errtheta = 0.01;   //   10mrad
	double errphi   = 0.01;   //   10mrad
	HepLorentzVector lvec = HepLorentzVector(pfo->getMomentum()[0], pfo->getMomentum()[1],pfo->getMomentum()[2],pfo->getEnergy());
	m_JetFitObject = new JetFitObject(lvec.e(), lvec.theta(), lvec.phi(),  1.0*std::sqrt(lvec.e()), errtheta, errphi);
}
//**********************************************
FSParticle::FSParticle(MCParticleImpl* pfo, string name, TLorentzVector p4):
		m_name(name), m_missed(false), m_fast(true)
{
	m_pfo                   = NULL;
	m_mcp                   = pfo;
	m_mass                  = Mass(m_name);
	m_rawFourMomentum       = p4;
	m_charge                = pfo->getCharge();
	m_Energy                = m_rawFourMomentum.E();
	m_pT                    = m_rawFourMomentum.Pt();
	m_pZ                    = m_rawFourMomentum.Pz();
	m_CosTheta              = m_rawFourMomentum.CosTheta();
	m_Rapidity              = m_rawFourMomentum.Rapidity();

	Int_t ntracks   = 0; 
	Int_t nclusters = 0;
	if ( fabs(m_charge)>0 ){
		ntracks++;
		m_trackId.push_back(pfo->id());
	}else{
		nclusters++;
		m_showerId.push_back(pfo->id());
	}
	double errtheta = 0.01;   //   10mrad
	double errphi   = 0.01;   //   10mrad
	HepLorentzVector lvec = HepLorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E());
	m_JetFitObject = new JetFitObject(lvec.e(), lvec.theta(), lvec.phi(), 
			1.0*std::sqrt(lvec.e()), errtheta, errphi);
}

//**********************************************
//**********************************************
FSParticle::~FSParticle()
{
	delete m_JetFitObject;
	FreeAll(m_trackId             );            
	FreeAll(m_showerId            );
	FreeAll(m_particleId          );
}
//**********************************************
//**********************************************
double FSParticle::sphericity() {

	// assumes I'm a flat jet (no vertex structure)
	TVector3 jetBoost = m_rawFourMomentum.BoostVector();
	TMatrixDSym sphMat(3);
	sphMat(0,0) = 0;
	sphMat(0,1) = 0;
	sphMat(0,2) = 0;
	sphMat(1,1) = 0;
	sphMat(1,2) = 0;
	sphMat(2,2) = 0;

	ReconstructedParticleVec tracks = m_pfo->getParticles();
	for (unsigned int i = 0; i<tracks.size(); i++) {
		ReconstructedParticle* v = tracks[i];
		TLorentzVector trkVec(v->getMomentum()[0],v->getMomentum()[1],v->getMomentum()[2],v->getEnergy());
		trkVec.Boost(-jetBoost);
		sphMat(0,0) += trkVec.X()*trkVec.X();
		sphMat(0,1) += trkVec.X()*trkVec.Y();
		sphMat(0,2) += trkVec.X()*trkVec.Z();
		sphMat(1,1) += trkVec.Y()*trkVec.Y();
		sphMat(1,2) += trkVec.Y()*trkVec.Z();
		sphMat(2,2) += trkVec.Z()*trkVec.Z();
	}
	double norm = sphMat(0,0)+sphMat(1,1)+sphMat(2,2);
	sphMat *= 1.5/norm;
	TVectorD eig;
	sphMat.EigenVectors(eig);
	//printf("eigen = %10.4f %10.4f %10.4f %10.4f\n",eig(0), eig(1), eig(2), eig(0)+eig(1)+eig(2));
	return eig(1)+eig(2);
}
//**********************************************
//**********************************************
int FSParticle::FindLeptonParent(MCParticle *mcp)
{
	int ret = 0;

	int PDGIDItself = abs(mcp->getPDG()); 
	if ( mcp->getParents().size() == 0 )      // original particle, no parent  
	{
		if ( PDGIDItself==11  ||  PDGIDItself==13 ){
		  	ret= PDGIDItself;
		}
	} 
	else if ( mcp->getParents().size() == 1 ) // Only one parent  
	{
		MCParticle *mcpa = mcp->getParents()[0];
		int PDGIDParent = abs(mcpa->getPDG()); 
		if ( PDGIDParent >100){
           PDGIDParent %= 1000;
		}
	   if ( ( PDGIDItself ==11 || PDGIDItself ==13 ) && (PDGIDParent == PDGIDItself) ){
			ret = FindLeptonParent(mcpa);	
		}else	
		if (     ( PDGIDItself ==11 || PDGIDItself ==13 ) 
				&& ( PDGIDParent ==13 || PDGIDParent ==15 || PDGIDParent ==23 || PDGIDParent ==24 || (PDGIDParent>100&&PDGIDParent<1000) )  
				){
			// parent is W/Z/tau/B/D and itself is a lepton 
			ret = PDGIDParent*100 + PDGIDItself;
		}
	}
	else                // More than one parent  
	{
	  printf("There are more than one parent: %2d \n", PDGIDItself);	
	}
	//
	return ret;
}


//**********************************************
//**********************************************
MCParticle * FSParticle::FindParton(MCParticle *mcp)
{
	MCParticle *ret = NULL;

	if ( mcp->getParents().size() == 0 )      // original particle, no parent  
	{
		if ( abs(mcp->getPDG())<=6 ){
		  	return mcp;
		}
		else return NULL;
	} 
	else if ( mcp->getParents().size() == 1 ) // Only one parent  
	{
		MCParticle *mcpa = mcp->getParents()[0];
		//
		if ( mcp->getPDG()==mcpa->getPDG() )
		{  
			// mother is same particle
			ret = FindParton(mcpa);	
		}
		else if ( abs(mcp->getPDG())<=6 && ( abs(mcpa->getPDG())>=23 && abs(mcpa->getPDG())<=25 ) )
		{
			// parent is W/Z and itself is a quark 
			ret = mcp;
		}
		else if ( mcp->getPDG()==21  && mcpa->getPDG()==25 )
		{
			// parent is Higgs and itself is a gluon 
			ret = mcp;
		}
		else if ( abs(mcp->getPDG())<= 6 && mcpa->getPDG()==21 )
		{
			// parent is gluon and itself is a quark, this kind of object will be neglected   
			ret = NULL;
		}	
		else if ( abs(mcpa->getPDG())<= 6 && mcp->getPDG()==21 )
		{
			// parent is quark and itself is a gluon 
			ret = FindParton(mcpa);	
		}	
		else if ( abs(mcpa->getPDG())>= 81 && abs(mcpa->getPDG())<=100 && mcpa->getParents().size()>=2 ) 
		{
			// parent is intermediate particle; need to check the MANY(>1) grandparents
			if ( abs(mcp->getPDG()) <= 6 || ( abs(mcp->getPDG()) >= 11 && abs(mcp->getPDG()) <= 16)  ) 
			{  
				// itself is a quark/lepton, jump over the intermediate paticle, search for the same quark/lepton, then ... 
				for ( unsigned int i=0; i<mcpa->getParents().size(); i++){
					MCParticle *mcpb = mcpa->getParents()[i];
					if ( mcpb->getPDG() == mcp->getPDG() ) {
						//printf( "here: %5d vs. %5d (%5d)\n", mcp->getPDG(),mcpb->getPDG(), mcpb->getParents()[0]->getPDG() );
						ret = FindParton(mcpb);
					}
				}
				if ( ret == NULL  && abs(mcp->getPDG()) <= 6 ) {
					for ( unsigned int i=0; i<mcpa->getParents().size(); i++ ) {
						MCParticle *mcpb = mcpa->getParents()[i];
						printf( "%5d vs. %5d, ", mcp->getPDG(),mcpb->getPDG() );
					}
					printf("\n");	
				}
			} 
			else 
			{  
				// itself is a hadron 
				// first; obtain partons from all grandparents
				map<MCParticle *, TVector3 > partonmap; partonmap.clear();
				map<MCParticle *, TVector3>::iterator it;
				for ( unsigned int i=0; i<mcpa->getParents().size(); i++){
					MCParticle *mcpb = mcpa->getParents()[i];
					MCParticle *mcpo = FindParton(mcpb);
					TVector3 v(mcpb->getMomentum());
					if ( mcpo != NULL) {
						it = partonmap.find(mcpo);
						if( it != partonmap.end() ) {
							v += it->second; // why add?
							//printf("the mother of %6d is %6d\n", (int)mcpb->getPDG(), (int)mcpo->getPDG());
						}
						partonmap[mcpo] = v;
					}
				}
				// then take the nearest one as the favorite parton
				double nearestangle = 2 * 3.14159;
				int selectedparton = -1, pdgid=0;
				MCParticle *pSelected = NULL;
				int j;
				TVector3 v2(mcp->getMomentum());
				for(j=0,it = partonmap.begin();it != partonmap.end(); it++,j++){
					if ( it->first ){
						double angle = v2.Angle(it->second);
						if(nearestangle > angle)
						{ 
							nearestangle = angle; selectedparton = j; pSelected = it->first; pdgid=(it->first)->getPDG();
						}
					}
				}
				ret = pSelected;
				if ( ret == NULL ) {
					//cout << partonmap.size() << " partons obtained." << endl;
					//cout << "Parton " << selectedparton << " selected. "<< nearestangle  << endl;
				}
				FreeAll(partonmap);
			}
		} 
		else
		{	
			ret = FindParton(mcpa);	
			if( 0 && ret==NULL){
				for (unsigned int k=0 ; k < mcpa->getParents().size() ; k++ ) {
					MCParticle 	* part    = mcpa->getParents()[k];
					printf("parent of %10d is %10d (%3d)\n", mcpa->getPDG(), part->getPDG(), (int)part->getParents().size() );
					for (unsigned int l=0 ; l < part->getParents().size() ; l++ ) {
						printf("parent %2d is %10d\n", l, part->getParents()[l]->getPDG() );
					}
				}
			}

		}	
	}
	else // More than one parent  
	{ 
		if ( abs(mcp->getPDG())<=6  ){
			MCParticle *mcpo2 = NULL;
			for( unsigned int i=0; i<mcp->getParents().size(); i++ ){
				MCParticle *mcpa = mcp->getParents()[i];
				if( mcp->getParents().size()==2 && abs(mcpa->getPDG()) == 11 ){
					mcpo2 = mcp;
					break;
				}else if( mcp->getPDG() == mcp->getPDG()  ){
					MCParticle *mcpo = FindParton(mcpa);
					mcpo2 = mcpo;
					break;
				}
			}
			ret=mcpo2;
		}
		else
		{	
			MCParticle *mcpo2 = NULL;
			for( unsigned int i=0; i<mcp->getParents().size(); i++ ){
				MCParticle *mcpa = mcp->getParents()[i];
				MCParticle *mcpo = FindParton(mcpa);
				if(mcpo2 && mcpo != mcpo2)
				{}
				//{cout << mcp->getPDG()<<" has Multiple parents other than intermidiates with different parton origin!" << endl;}
				mcpo2 = mcpo;
			}
			//
			ret = mcpo2;
		}
	}
	//
	return ret;
}
//**********************************************
//
//   UTILITY FUNCTIONS
//
//**********************************************
TLorentzVector
FSParticle::fourMomentum( ReconstructedParticle* sh ){
	return TLorentzVector(sh->getMomentum()[0],sh->getMomentum()[1],sh->getMomentum()[2],sh->getEnergy());
}
TLorentzVector
FSParticle::fourMomentum( MCParticle* sh ){
	return TLorentzVector(sh->getMomentum()[0],sh->getMomentum()[1],sh->getMomentum()[2],sh->getEnergy());
}

bool
FSParticle::duplicate(FSParticle* fsp, int full ){
	/*
	//
	vector<int> trackId2  = fsp->trackId();
	vector<int> showerId2 = fsp->showerId();
	//
	if (m_trackId.size() > 0 && trackId2.size() > 0){
		for (unsigned int i = 0; i < m_trackId.size(); i++){
			for (unsigned int j = 0; j < trackId2.size(); j++){
				if (m_trackId[i] == trackId2[j]) return true;
			}
		}
		if ( name() == fsp->name() && m_trackId[0] > trackId2[0] ) return true;
	}
	// if there are two jets, one has only tracks and the other only showers, how to check duplication by track and shower information? 
	if (m_showerId.size() > 0 && showerId2.size() > 0){
		for (unsigned int i = 0; i < m_showerId.size(); i++){
			for (unsigned int j = 0; j < showerId2.size(); j++){
				if (m_showerId[i] == showerId2[j]) return true;
			}
		}
		//if ( name()=="jet"   && name() == fsp->name() && m_showerId[0] > showerId2[0]) return true;
		if ( name()=="gamma" && name() == fsp->name() && m_showerId[0] > showerId2[0]) return true;
	}
	*/
	//
	if( full==1 ){	
		/*
		vector<int> particleId2 = fsp->particleId();
		if (m_particleId.size() > 0 && particleId2.size() > 0){
			for (unsigned int i = 0; i < m_particleId.size(); i++){
				for (unsigned int j = 0; j < particleId2.size(); j++){
					if (m_particleId[i] == particleId2[j]) return true;
				}
			}
			if ( name() == fsp->name() && m_particleId[0] > particleId2[0]) return true;
		}
		*/
		vector<int> trackId2  = fsp->trackId();
		vector<int> showerId2 = fsp->showerId();
		//
		if (m_trackId.size() > 0 && trackId2.size() > 0){
			for (unsigned int i = 0; i < m_trackId.size(); i++){
				for (unsigned int j = 0; j < trackId2.size(); j++){
					if (m_trackId[i] == trackId2[j]) return true;
				}
			}
			if ( name() == fsp->name() && m_trackId[0] > trackId2[0] ) return true;
		}
		// if there are two jets, one has only tracks and the other only showers, how to check duplication by track and shower information? 
		if (m_showerId.size() > 0 && showerId2.size() > 0){
			for (unsigned int i = 0; i < m_showerId.size(); i++){
				for (unsigned int j = 0; j < showerId2.size(); j++){
					if (m_showerId[i] == showerId2[j]) return true;
				}
			}
			//if ( name()=="jet"   && name() == fsp->name() && m_showerId[0] > showerId2[0]) return true;
			if ( name()=="gamma" && name() == fsp->name() && m_showerId[0] > showerId2[0]) return true;
		}
	}else{
		//
		vector<int> particleId2 = fsp->particleId();
		//cout<<"duplicate "<<m_particleId.size() << particleId2.size() <<endl;; 
		//
		if (m_particleId.size() > 0 && particleId2.size() > 0){
			for (unsigned int i = 0; i < m_particleId.size(); i++){
				for (unsigned int j = 0; j < particleId2.size(); j++){
					if (m_particleId[i] == particleId2[j]) return true;
				}
			}
			if ( name() == fsp->name() && m_particleId[0] > particleId2[0]) return true;
		}
		/*
		   vector<int> trackId2 = fsp->trackId();
		if (m_trackId.size() > 0 && trackId2.size() > 0){
			for (unsigned int i = 0; i < m_trackId.size(); i++){
				for (unsigned int j = 0; j < trackId2.size(); j++){
					if (m_trackId[i] == trackId2[j]) return true;
				}
			}
			if ( name() == fsp->name() && m_trackId[0] > trackId2[0]) return true;
		}

		vector<int> showerId2 = fsp->showerId();
		if (m_showerId.size() > 0 && showerId2.size() > 0){
			for (unsigned int i = 0; i < m_showerId.size(); i++){
				for (unsigned int j = 0; j < showerId2.size(); j++){
					if (m_showerId[i] == showerId2[j]) return true;
				}
			}
			if ( name()=="gamma" && name() == fsp->name() && m_showerId[0] > showerId2[0]) return true;
		}
		*/
	}
	//
	//
	return false;
}

void 
FSParticle::PrintTrackAndShowers(){
	cout<<" Tracks: "; 
	for(unsigned int i=0; i<m_trackId .size(); i++) cout<<m_trackId [i]<<", "; 
	cout<<";Showers "; 
	for(unsigned int i=0; i<m_showerId.size(); i++) cout<<m_showerId[i]<<", "; 
	cout<<endl; 
}

//**********************************************
//
//   PARSE FINAL STATE INFORMATION
//
//**********************************************
FSInfo::FSInfo(string FSName, NTupleHelper* nt, NTupleHelper* ntgen){
	m_FSName            = FSName;
	m_NT                = nt;
	m_NTGen             = ntgen;
	m_nMissingParticles = 0;
	m_missingMassFit    = 0;
	m_missingMassValue  =-1; 
	m_Constrain4Mom     = 1;
	m_FSCuts.clear();

	cout << "\t\t"<< endl;
	cout << "FSClasser:  Initializing Final State " << FSName << endl;

	m_missedParticle = "";
	string singleParticles[8] = {
		"e-"     ,
		"e+"     ,
		"mu-"    ,
		"mu+"    ,
		"tau-"   ,
		"tau+"   ,
		"gamma"  ,
		//
		 "jet"    
	};
	//
	int index = 0, idash=0;
	vector<string> particleNamesTmp;
	vector<int>    particleStatusTmp;
	for (int i = m_FSName.size()-1; i >= 0 && index < 11; i--){
		string digit = m_FSName.substr(i,1);
		if (            (digit == "0") || (digit == "1") || 
				(digit == "2") || (digit == "3") || 
				(digit == "4") || (digit == "5") || 
				(digit == "6") || (digit == "7") || 
				(digit == "8") || (digit == "9")
				){
			int num = atoi(digit.c_str());
			for (int j = 0; j < num; j++){
				particleNamesTmp .push_back(singleParticles[index]);
				particleStatusTmp.push_back(1);
			}
			index++;
		}else if (digit == "M" || digit == "m" ){
			m_nMissingParticles++;
			m_missingMassFit = true; 
			particleNamesTmp .push_back(singleParticles[index]);
			particleStatusTmp.push_back(0);
			m_missedParticle = singleParticles[index];
			setMissingMass(Mass(m_missedParticle));
			index++;
		}else if (digit == "_" ){
			if( idash==0){
				index =  7;
				idash++;
			}
		}
	}
	//
	for (int i = particleNamesTmp.size()-1; i >= 0; i--){
		m_particleNames .push_back(particleNamesTmp[i] );
		m_particleStatus.push_back(particleStatusTmp[i]);
	}
	//
	m_nChargedParticles = 0;
	m_decayCode1 = 0;
	m_decayCode2 = 0;
	for (unsigned int i = 0; i < m_particleNames.size(); i++){
		if      (m_particleNames[i] == "e-"      ){ m_decayCode1 += 1;         if(m_particleStatus[i]!=0)m_nChargedParticles += 1; }
		else if (m_particleNames[i] == "e+"      ){ m_decayCode1 += 10;        if(m_particleStatus[i]!=0)m_nChargedParticles += 2; }
		else if (m_particleNames[i] == "mu-"     ){ m_decayCode1 += 100;       if(m_particleStatus[i]!=0)m_nChargedParticles += 1; }
		else if (m_particleNames[i] == "mu+"     ){ m_decayCode1 += 1000;      if(m_particleStatus[i]!=0)m_nChargedParticles += 1; }
		else if (m_particleNames[i] == "tau-"    ){ m_decayCode1 += 10000;     if(m_particleStatus[i]!=0)m_nChargedParticles += 0; }
		else if (m_particleNames[i] == "tau+"    ){ m_decayCode1 += 100000;    if(m_particleStatus[i]!=0)m_nChargedParticles += 0; }
		else if (m_particleNames[i] == "gamma"   ){ m_decayCode1 += 1000000;   if(m_particleStatus[i]!=0)m_nChargedParticles += 1; }
		else if (m_particleNames[i] ==  "jet"    ){ m_decayCode2 += 1;         if(m_particleStatus[i]!=0)m_nChargedParticles += 0; }
	}
	//
	if( inclusive() )          m_Constrain4Mom=0;
	if(m_nMissingParticles>1)  m_Constrain4Mom=0;
	//
	Print();
	/*
	for (unsigned int i = 0; i < m_particleNames.size(); i++){
		cout << "FSClasser:      " << m_particleNames[i];
		if( m_particleStatus[i] ) cout << endl;
		else                      cout <<" is missing"<<endl;
	}
	*/
	FreeAll(particleNamesTmp);
	FreeAll(particleStatusTmp);
}

FSInfo::~FSInfo(){
	FreeAll   (  m_particleNames );
	FreeAll   (  m_particleStatus);
	FreeDelAll(  m_FSCuts        );
	if(m_NT    ) delete m_NT   ;
	if(m_NTGen ) delete m_NTGen;
}

//**********************************************
//
//   FSCut constructor
//
//**********************************************

FSCut::FSCut(const vector<string> initialization){


	// divide the initialization string into seperate words

	//vector<string> words = FSInfo::parseString(initialization);
	vector<string> words = initialization; 

	if (words.size() != 5){
		cout << "FSClasser ERROR: wrong arguments to FSCut parameter: " << endl;
		cout << initialization[0] << endl;
		exit(0);
	}

	// save input as member data

	m_FSName       = words[0];
	m_submodeName  = words[1];
	m_cutType      = words[2];
	m_lowCut  = atof(words[3].c_str());
	m_highCut = atof(words[4].c_str());


	// quick checks

	if (            m_cutType != "RawRecoil"        && m_cutType != "RawMass"        &&
			m_cutType != "FitRecoil"        && m_cutType != "FitMass"        &&
			m_cutType != "RawRecoilSquared" && m_cutType != "RawMassSquared" &&
			m_cutType != "FitRecoilSquared" && m_cutType != "FitMassSquared" ){
		cout << "FSClasser ERROR: wrong arguments to FSCut parameter: " << endl;
		cout << initialization[0] << endl;
		cout << "cutType must be RawRecoil, RawMass, IntRecoil, FitRecoil, etc. " << endl;
		exit(0);
	}

}

//**********************************************
//
//   FSInfo functions to evaluate FSCuts
//
//**********************************************
vector< vector<unsigned int> >&
FSInfo::submodeIndices(const string& submodeName){

	if (m_submodeIndices.find(submodeName) != m_submodeIndices.end())
		return m_submodeIndices[submodeName];

	vector<string> submodeParticles = FSInfo::getParticleNamesFromFSName( submodeName );

	static vector< vector<unsigned int> > indices;

	for (unsigned int i = 0; i < submodeParticles.size(); i++){

		vector< vector<unsigned int> > indicesTemp = indices;
		indices.clear();

		vector<unsigned int> pList;
		for (unsigned int j = 0; j < m_particleNames.size(); j++){
			if (submodeParticles[i] == m_particleNames[j]) pList.push_back(j);
		}
		if (pList.size() == 0) return indices;

		for (unsigned int ipl = 0; ipl < pList.size(); ipl++){
			if (i == 0){
				vector<unsigned int> combo;
				combo.push_back(pList[ipl]);
				indices.push_back(combo);
			}
			else{
				for (unsigned int itmp = 0; itmp < indicesTemp.size(); itmp++){
					vector<unsigned int> combo = indicesTemp[itmp];
					bool duplicate = false;
					for (unsigned int ic = 0; ic < combo.size(); ic++){
						if (pList[ipl] <= combo[ic]){
							duplicate = true;
							continue;
						}
					}
					if (!duplicate){
						combo.push_back(pList[ipl]);
						indices.push_back(combo);
					}
				}
			}
		}
	}

	m_submodeIndices[submodeName] = indices;
	return m_submodeIndices[submodeName];

}

TVector3 
FSInfo::HiggsBoostVector( const vector<FSParticle*>& particleCombination, string RecoilSide ){

	TVector3 v3(0,0,0); 
	vector< vector<unsigned int> > indices = submodeIndices(RecoilSide);
	if( indices.size()>0){
		TLorentzVector v4(0.0, 0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < indices.size(); i++){
			vector<unsigned int> indexCombo = indices[i];

			TLorentzVector p4(0.0, 0.0, 0.0, 0.0);
			for (unsigned int j = 0; j < indexCombo.size(); j++){
				p4 += particleCombination[indexCombo[j]]->rawFourMomentum();
			}
			if ( p4.Rho()>v4.Rho()) v4=p4; 

		}

		double Beta   = v4.Beta();
		v3 = v4.Vect();
		v3.SetMag(Beta);
	}
	return v3;

}



bool
FSInfo::evaluateFSCuts( const vector<FSParticle*>& particleCombination,
		const TLorentzVector& pInitial, string fourVectorType ){



	for (unsigned int icut = 0; icut < m_FSCuts.size(); icut++){
		FSCut* fscut = m_FSCuts[icut];
		bool lRaw = fscut->Raw();
		bool lFit = fscut->Fit();
		if ((lRaw && fourVectorType != "Raw") ||
				(lFit && fourVectorType != "Fit")) continue;
		bool pass = false;

		vector< vector<unsigned int> > indices = submodeIndices(fscut->submodeName());

		for (unsigned int i = 0; i < indices.size(); i++){
			vector<unsigned int> indexCombo = indices[i];

			TLorentzVector pTotal(0.0, 0.0, 0.0, 0.0);
			for (unsigned int j = 0; j < indexCombo.size(); j++){
				if (lRaw) pTotal += particleCombination[indexCombo[j]]->rawFourMomentum();
				if (lFit) pTotal += particleCombination[indexCombo[j]]->fitFourMomentum();
			}

			if (fscut->Mass()){ 
				double x;
				if (fscut->Squared()){ x = pTotal.M2(); }
				else                 { x = pTotal.M(); }
				if (x > fscut->lowCut() && x < fscut->highCut()){ pass = true; break; }
			}

			else if (fscut->Recoil()){ 
				TLorentzVector pMiss = pInitial - pTotal;
				double x;
				if (fscut->Squared()){ x = pMiss.M2(); }
				else                 { x = pMiss.M(); }
				if (x > fscut->lowCut() && x < fscut->highCut()){ pass = true; break; }
			}
		}


		if (!pass) return false;

	}

	return true;

}


//**********************************************
//
//   FSInfo: MODE NUMBERING UTILITIES, ETC.
//
//**********************************************
vector<string>
FSInfo::getParticleNamesFromFSName( const string& FSName ){


	// some quick checks

	if ((FSName.size() == 0) ||
			(FSName.find("_") == string::npos)){
		cout << "FSClasser ERROR: error in final state name: " << FSName << endl;
		exit(1);
	}


	// a list of allowed particle names

	string singleParticles[8] = {
		"e-"   ,  
		"e+"   ,
		"mu-"  , 
		"mu+"  ,
		"tau-" , 
		"tau+" ,
		"gamma", 
		"jet" 
	};

	// parse FSName digit by digit, starting at the end

	int index = 0, idash=0;
	vector<string> particleNamesTmp;
	for (int i = FSName.size()-1; i >= 0 && index < 11; i--){
		string digit = FSName.substr(i,1);
		if (  (digit == "0") || (digit == "1") || 
				(digit == "2") || (digit == "3") || 
				(digit == "4") || (digit == "5") || 
				(digit == "6") || (digit == "7") || 
				(digit == "8") || (digit == "9")){
			int num = atoi(digit.c_str());
			for (int j = 0; j < num; j++){
				particleNamesTmp.push_back(singleParticles[index]);
			}
			index++;
		}else if (digit == "M" || digit == "m" ){
			particleNamesTmp .push_back(singleParticles[index]);
			index++;
		}else if (digit == "_" ){
			if( idash==0){
				index = 7;
				idash++;
			}
		}
		else{
			break;
		}
	}

	// make sure we have particles

	if (particleNamesTmp.size() == 0){
		cout << "FSClasser ERROR: error in final state name: " << FSName << endl;
		exit(1);
	}

	// now reverse the order of the particle names

	vector<string> particleNames;
	for (int i = particleNamesTmp.size()-1; i >= 0; i--){
		particleNames   .push_back(particleNamesTmp[i]);
	}

	return particleNames;

}

vector<int>
FSInfo::getDecayCodeFromParticleNames( const vector<string>& particleNames ){

	int decayCode1 = 0;
	int decayCode2 = 0;
	for (unsigned int i = 0; i < particleNames.size(); i++){
		if      (particleNames[i] == "e-"       ){ decayCode1 += 1;           }
		else if (particleNames[i] == "e+"       ){ decayCode1 += 10;          }
		else if (particleNames[i] == "mu-"      ){ decayCode1 += 1000;        }
		else if (particleNames[i] == "mu+"      ){ decayCode1 += 10000;       }
		else if (particleNames[i] == "tau-"     ){ decayCode1 += 100000;      }
		else if (particleNames[i] == "tau+"     ){ decayCode1 += 1000000;     }
		else if (particleNames[i] == "gamma"    ){ decayCode2 += 10000000;    }
		//
		else if (particleNames[i] ==  "jet"     ){ decayCode2 += 1;           }
	}
	vector<int> decayCode; decayCode.clear();
	decayCode.push_back(decayCode1);
	decayCode.push_back(decayCode2);
	return decayCode;
}


vector<int>
FSInfo::getDecayCodeFromFSName( const string& FSName ){

	return getDecayCodeFromParticleNames(getParticleNamesFromFSName(FSName));

}


int
FSInfo::getNChargedParticlesFromParticleNames( const vector<string>& particleNames, const vector<string>& particleStatus  ){

	int nChargedParticles = 0;
	for (unsigned int i = 0; i < particleNames.size(); i++){
		if      (particleNames[i] ==  "jet"   ){ if(particleStatus[i]!=0) nChargedParticles += 0; }
		else if (particleNames[i] == "gamma"  ){ if(particleStatus[i]!=0) nChargedParticles += 0; }
		else if (particleNames[i] == "e-"     ){ if(particleStatus[i]!=0) nChargedParticles += 1; }
		else if (particleNames[i] == "e+"     ){ if(particleStatus[i]!=0) nChargedParticles += 1; }
		else if (particleNames[i] == "mu-"    ){ if(particleStatus[i]!=0) nChargedParticles += 1; }
		else if (particleNames[i] == "mu+"    ){ if(particleStatus[i]!=0) nChargedParticles += 1; }
		else if (particleNames[i] == "tau-"   ){ if(particleStatus[i]!=0) nChargedParticles += 1; }
		else if (particleNames[i] == "tau+"   ){ if(particleStatus[i]!=0) nChargedParticles += 1; }
	}
	return nChargedParticles;

}

vector<string>
FSInfo::parseString( const string& inputString ){

	vector<string> words;
	string word("");
	for (unsigned int j = 0; j < inputString.size(); j++){
		if (!isspace(inputString[j])){
			word += inputString[j];
			if ((j == (inputString.size()-1))&&(!word.empty())){
				words.push_back(word);
				word = "";
			}
		} 
		else if (!word.empty()){
			words.push_back(word);
			word = "";
		}
	}
	return words;
}
