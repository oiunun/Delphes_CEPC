//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Nov 19 22:08:48 2022 by ROOT version 6.22/06
// from TTree Delphes/Analysis tree
// found on file: ../../rootfile/B.root
//////////////////////////////////////////////////////////

#ifndef b_h
#define b_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

class b {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxEvent = 1;
   static constexpr Int_t kMaxWeight = 1;
   static constexpr Int_t kMaxParticle = 282;
   static constexpr Int_t kMaxGenVertex = 14;
   static constexpr Int_t kMaxTrack = 43;
   static constexpr Int_t kMaxEFlowTrack = 40;
   static constexpr Int_t kMaxParticleFlowCandidate = 60;
   static constexpr Int_t kMaxElectron = 2;
   static constexpr Int_t kMaxPhoton = 9;
   static constexpr Int_t kMaxMuon = 7;

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   Long64_t        Event_Number[kMaxEvent];   //[Event_]
   Float_t         Event_ReadTime[kMaxEvent];   //[Event_]
   Float_t         Event_ProcTime[kMaxEvent];   //[Event_]
   Int_t           Event_ProcessID[kMaxEvent];   //[Event_]
   Int_t           Event_MPI[kMaxEvent];   //[Event_]
   Float_t         Event_Weight[kMaxEvent];   //[Event_]
   Float_t         Event_CrossSection[kMaxEvent];   //[Event_]
   Float_t         Event_CrossSectionError[kMaxEvent];   //[Event_]
   Float_t         Event_Scale[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQED[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQCD[kMaxEvent];   //[Event_]
   Int_t           Event_ID1[kMaxEvent];   //[Event_]
   Int_t           Event_ID2[kMaxEvent];   //[Event_]
   Float_t         Event_X1[kMaxEvent];   //[Event_]
   Float_t         Event_X2[kMaxEvent];   //[Event_]
   Float_t         Event_ScalePDF[kMaxEvent];   //[Event_]
   Float_t         Event_PDF1[kMaxEvent];   //[Event_]
   Float_t         Event_PDF2[kMaxEvent];   //[Event_]
   Int_t           Event_size;
   Int_t           Weight_;
   UInt_t          Weight_fUniqueID[kMaxWeight];   //[Weight_]
   UInt_t          Weight_fBits[kMaxWeight];   //[Weight_]
   Float_t         Weight_Weight[kMaxWeight];   //[Weight_]
   Int_t           Weight_size;
   Int_t           Particle_;
   UInt_t          Particle_fUniqueID[kMaxParticle];   //[Particle_]
   UInt_t          Particle_fBits[kMaxParticle];   //[Particle_]
   Int_t           Particle_PID[kMaxParticle];   //[Particle_]
   Int_t           Particle_Status[kMaxParticle];   //[Particle_]
   Int_t           Particle_IsPU[kMaxParticle];   //[Particle_]
   Int_t           Particle_M1[kMaxParticle];   //[Particle_]
   Int_t           Particle_M2[kMaxParticle];   //[Particle_]
   Int_t           Particle_D1[kMaxParticle];   //[Particle_]
   Int_t           Particle_D2[kMaxParticle];   //[Particle_]
   Int_t           Particle_Charge[kMaxParticle];   //[Particle_]
   Float_t         Particle_Mass[kMaxParticle];   //[Particle_]
   Float_t         Particle_E[kMaxParticle];   //[Particle_]
   Float_t         Particle_Px[kMaxParticle];   //[Particle_]
   Float_t         Particle_Py[kMaxParticle];   //[Particle_]
   Float_t         Particle_Pz[kMaxParticle];   //[Particle_]
   Float_t         Particle_P[kMaxParticle];   //[Particle_]
   Float_t         Particle_PT[kMaxParticle];   //[Particle_]
   Float_t         Particle_Eta[kMaxParticle];   //[Particle_]
   Float_t         Particle_Phi[kMaxParticle];   //[Particle_]
   Float_t         Particle_Rapidity[kMaxParticle];   //[Particle_]
   Float_t         Particle_T[kMaxParticle];   //[Particle_]
   Float_t         Particle_X[kMaxParticle];   //[Particle_]
   Float_t         Particle_Y[kMaxParticle];   //[Particle_]
   Float_t         Particle_Z[kMaxParticle];   //[Particle_]
   Int_t           Particle_size;
   Int_t           GenVertex_;
   UInt_t          GenVertex_fUniqueID[kMaxGenVertex];   //[GenVertex_]
   UInt_t          GenVertex_fBits[kMaxGenVertex];   //[GenVertex_]
   Float_t         GenVertex_T[kMaxGenVertex];   //[GenVertex_]
   Float_t         GenVertex_X[kMaxGenVertex];   //[GenVertex_]
   Float_t         GenVertex_Y[kMaxGenVertex];   //[GenVertex_]
   Float_t         GenVertex_Z[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_ErrorT[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_ErrorX[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_ErrorY[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_ErrorZ[kMaxGenVertex];   //[GenVertex_]
   Int_t           GenVertex_Index[kMaxGenVertex];   //[GenVertex_]
   Int_t           GenVertex_NDF[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_Sigma[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_SumPT2[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_GenSumPT2[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_GenDeltaZ[kMaxGenVertex];   //[GenVertex_]
   Double_t        GenVertex_BTVSumPT2[kMaxGenVertex];   //[GenVertex_]
   TRefArray       GenVertex_Constituents[kMaxGenVertex];
   Int_t           GenVertex_size;
   Int_t           Track_;
   UInt_t          Track_fUniqueID[kMaxTrack];   //[Track_]
   UInt_t          Track_fBits[kMaxTrack];   //[Track_]
   Int_t           Track_PID[kMaxTrack];   //[Track_]
   Double_t        Track_Chi_k[kMaxTrack];   //[Track_]
   Double_t        Track_Chi_pi[kMaxTrack];   //[Track_]
   Double_t        Track_Counting_eff[kMaxTrack];   //[Track_]
   Double_t        Track_L_DC[kMaxTrack];   //[Track_]
   Int_t           Track_Truth_PID[kMaxTrack];   //[Track_]
   Float_t         Track_Truth_P[kMaxTrack];   //[Track_]
   Float_t         Track_Truth_CosTheta[kMaxTrack];   //[Track_]
   Float_t         Track_Truth_Phi[kMaxTrack];   //[Track_]
   Int_t           Track_Charge[kMaxTrack];   //[Track_]
   Float_t         Track_P[kMaxTrack];   //[Track_]
   Float_t         Track_PT[kMaxTrack];   //[Track_]
   Float_t         Track_Eta[kMaxTrack];   //[Track_]
   Float_t         Track_Phi[kMaxTrack];   //[Track_]
   Float_t         Track_CtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_CosTheta[kMaxTrack];   //[Track_]
   Float_t         Track_C[kMaxTrack];   //[Track_]
   Float_t         Track_Mass[kMaxTrack];   //[Track_]
   Float_t         Track_Prob_K[kMaxTrack];   //[Track_]
   Float_t         Track_Prob_Pi[kMaxTrack];   //[Track_]
   Float_t         Track_Prob_P[kMaxTrack];   //[Track_]
   Float_t         Track_EtaOuter[kMaxTrack];   //[Track_]
   Float_t         Track_PhiOuter[kMaxTrack];   //[Track_]
   Float_t         Track_T[kMaxTrack];   //[Track_]
   Float_t         Track_X[kMaxTrack];   //[Track_]
   Float_t         Track_Y[kMaxTrack];   //[Track_]
   Float_t         Track_Z[kMaxTrack];   //[Track_]
   Float_t         Track_TOuter[kMaxTrack];   //[Track_]
   Float_t         Track_XOuter[kMaxTrack];   //[Track_]
   Float_t         Track_YOuter[kMaxTrack];   //[Track_]
   Float_t         Track_ZOuter[kMaxTrack];   //[Track_]
   Float_t         Track_Xd[kMaxTrack];   //[Track_]
   Float_t         Track_Yd[kMaxTrack];   //[Track_]
   Float_t         Track_Zd[kMaxTrack];   //[Track_]
   Float_t         Track_L[kMaxTrack];   //[Track_]
   Float_t         Track_D0[kMaxTrack];   //[Track_]
   Float_t         Track_DZ[kMaxTrack];   //[Track_]
   Float_t         Track_Nclusters[kMaxTrack];   //[Track_]
   Float_t         Track_Nclusters_err[kMaxTrack];   //[Track_]
   Float_t         Track_dNdx[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorP[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPT[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhi[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorCtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorT[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorDZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorC[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0Phi[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0C[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0DZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0CtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhiC[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhiDZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhiCtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorCDZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorCCtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorDZCtgTheta[kMaxTrack];   //[Track_]
   TRef            Track_Particle[kMaxTrack];
   Int_t           Track_VertexIndex[kMaxTrack];   //[Track_]
   Int_t           Track_size;
   Int_t           EFlowTrack_;
   UInt_t          EFlowTrack_fUniqueID[kMaxEFlowTrack];   //[EFlowTrack_]
   UInt_t          EFlowTrack_fBits[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_PID[kMaxEFlowTrack];   //[EFlowTrack_]
   Double_t        EFlowTrack_Chi_k[kMaxEFlowTrack];   //[EFlowTrack_]
   Double_t        EFlowTrack_Chi_pi[kMaxEFlowTrack];   //[EFlowTrack_]
   Double_t        EFlowTrack_Counting_eff[kMaxEFlowTrack];   //[EFlowTrack_]
   Double_t        EFlowTrack_L_DC[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_Truth_PID[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Truth_P[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Truth_CosTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Truth_Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_Charge[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_P[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_PT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Eta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_CtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_CosTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_C[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Mass[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Prob_K[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Prob_Pi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Prob_P[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_EtaOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_PhiOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_T[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_X[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Y[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Z[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_TOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_XOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_YOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ZOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Xd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Yd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Zd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_L[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_D0[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_DZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Nclusters[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Nclusters_err[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_dNdx[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorP[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorC[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0C[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0DZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0CtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiC[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorDZCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   TRef            EFlowTrack_Particle[kMaxEFlowTrack];
   Int_t           EFlowTrack_VertexIndex[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_size;
   Int_t           ParticleFlowCandidate_;
   UInt_t          ParticleFlowCandidate_fUniqueID[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   UInt_t          ParticleFlowCandidate_fBits[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Int_t           ParticleFlowCandidate_PID[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Int_t           ParticleFlowCandidate_Charge[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_E[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_P[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_PT[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Eta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Phi[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_CtgTheta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_C[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Mass[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_EtaOuter[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_PhiOuter[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_T[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_X[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Y[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Z[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_TOuter[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_XOuter[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_YOuter[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ZOuter[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Xd[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Yd[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Zd[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_L[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_D0[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_DZ[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Nclusters[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Nclusters_err[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_dNdx[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorP[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorPT[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorPhi[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorCtgTheta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorT[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorD0[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorDZ[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorC[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorD0Phi[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorD0C[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorD0DZ[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorD0CtgTheta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorPhiC[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorPhiDZ[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorPhiCtgTheta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorCDZ[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorCCtgTheta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_ErrorDZCtgTheta[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Int_t           ParticleFlowCandidate_VertexIndex[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Int_t           ParticleFlowCandidate_NTimeHits[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Eem[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Ehad[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Etrk[kMaxParticleFlowCandidate];   //[ParticleFlowCandidate_]
   Float_t         ParticleFlowCandidate_Edges[kMaxParticleFlowCandidate][4];   //[ParticleFlowCandidate_]
   TRefArray       ParticleFlowCandidate_Particles[kMaxParticleFlowCandidate];
   Int_t           ParticleFlowCandidate_size;
   Int_t           Electron_;
   UInt_t          Electron_fUniqueID[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fBits[kMaxElectron];   //[Electron_]
   Float_t         Electron_PT[kMaxElectron];   //[Electron_]
   Float_t         Electron_Eta[kMaxElectron];   //[Electron_]
   Float_t         Electron_Phi[kMaxElectron];   //[Electron_]
   Float_t         Electron_T[kMaxElectron];   //[Electron_]
   Int_t           Electron_Charge[kMaxElectron];   //[Electron_]
   Float_t         Electron_EhadOverEem[kMaxElectron];   //[Electron_]
   TRef            Electron_Particle[kMaxElectron];
   Float_t         Electron_IsolationVar[kMaxElectron];   //[Electron_]
   Float_t         Electron_IsolationVarRhoCorr[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtCharged[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtNeutral[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtChargedPU[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPt[kMaxElectron];   //[Electron_]
   Float_t         Electron_D0[kMaxElectron];   //[Electron_]
   Float_t         Electron_DZ[kMaxElectron];   //[Electron_]
   Float_t         Electron_ErrorD0[kMaxElectron];   //[Electron_]
   Float_t         Electron_ErrorDZ[kMaxElectron];   //[Electron_]
   Int_t           Electron_size;
   Int_t           Photon_;
   UInt_t          Photon_fUniqueID[kMaxPhoton];   //[Photon_]
   UInt_t          Photon_fBits[kMaxPhoton];   //[Photon_]
   Float_t         Photon_PT[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Eta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Phi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_E[kMaxPhoton];   //[Photon_]
   Float_t         Photon_T[kMaxPhoton];   //[Photon_]
   Float_t         Photon_EhadOverEem[kMaxPhoton];   //[Photon_]
   TRefArray       Photon_Particles[kMaxPhoton];
   Float_t         Photon_IsolationVar[kMaxPhoton];   //[Photon_]
   Float_t         Photon_IsolationVarRhoCorr[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtCharged[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtNeutral[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtChargedPU[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPt[kMaxPhoton];   //[Photon_]
   Float_t         Photon_CosTheta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Truth_Phi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Truth_E[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Truth_CosTheta[kMaxPhoton];   //[Photon_]
   Int_t           Photon_Status[kMaxPhoton];   //[Photon_]
   Int_t           Photon_size;
   Int_t           Muon_;
   UInt_t          Muon_fUniqueID[kMaxMuon];   //[Muon_]
   UInt_t          Muon_fBits[kMaxMuon];   //[Muon_]
   Float_t         Muon_PT[kMaxMuon];   //[Muon_]
   Float_t         Muon_Eta[kMaxMuon];   //[Muon_]
   Float_t         Muon_Phi[kMaxMuon];   //[Muon_]
   Float_t         Muon_T[kMaxMuon];   //[Muon_]
   Int_t           Muon_Charge[kMaxMuon];   //[Muon_]
   TRef            Muon_Particle[kMaxMuon];
   Float_t         Muon_IsolationVar[kMaxMuon];   //[Muon_]
   Float_t         Muon_IsolationVarRhoCorr[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPtCharged[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPtNeutral[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPtChargedPU[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPt[kMaxMuon];   //[Muon_]
   Float_t         Muon_D0[kMaxMuon];   //[Muon_]
   Float_t         Muon_DZ[kMaxMuon];   //[Muon_]
   Float_t         Muon_ErrorD0[kMaxMuon];   //[Muon_]
   Float_t         Muon_ErrorDZ[kMaxMuon];   //[Muon_]
   Int_t           Muon_size;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_ReadTime;   //!
   TBranch        *b_Event_ProcTime;   //!
   TBranch        *b_Event_ProcessID;   //!
   TBranch        *b_Event_MPI;   //!
   TBranch        *b_Event_Weight;   //!
   TBranch        *b_Event_CrossSection;   //!
   TBranch        *b_Event_CrossSectionError;   //!
   TBranch        *b_Event_Scale;   //!
   TBranch        *b_Event_AlphaQED;   //!
   TBranch        *b_Event_AlphaQCD;   //!
   TBranch        *b_Event_ID1;   //!
   TBranch        *b_Event_ID2;   //!
   TBranch        *b_Event_X1;   //!
   TBranch        *b_Event_X2;   //!
   TBranch        *b_Event_ScalePDF;   //!
   TBranch        *b_Event_PDF1;   //!
   TBranch        *b_Event_PDF2;   //!
   TBranch        *b_Event_size;   //!
   TBranch        *b_Weight_;   //!
   TBranch        *b_Weight_fUniqueID;   //!
   TBranch        *b_Weight_fBits;   //!
   TBranch        *b_Weight_Weight;   //!
   TBranch        *b_Weight_size;   //!
   TBranch        *b_Particle_;   //!
   TBranch        *b_Particle_fUniqueID;   //!
   TBranch        *b_Particle_fBits;   //!
   TBranch        *b_Particle_PID;   //!
   TBranch        *b_Particle_Status;   //!
   TBranch        *b_Particle_IsPU;   //!
   TBranch        *b_Particle_M1;   //!
   TBranch        *b_Particle_M2;   //!
   TBranch        *b_Particle_D1;   //!
   TBranch        *b_Particle_D2;   //!
   TBranch        *b_Particle_Charge;   //!
   TBranch        *b_Particle_Mass;   //!
   TBranch        *b_Particle_E;   //!
   TBranch        *b_Particle_Px;   //!
   TBranch        *b_Particle_Py;   //!
   TBranch        *b_Particle_Pz;   //!
   TBranch        *b_Particle_P;   //!
   TBranch        *b_Particle_PT;   //!
   TBranch        *b_Particle_Eta;   //!
   TBranch        *b_Particle_Phi;   //!
   TBranch        *b_Particle_Rapidity;   //!
   TBranch        *b_Particle_T;   //!
   TBranch        *b_Particle_X;   //!
   TBranch        *b_Particle_Y;   //!
   TBranch        *b_Particle_Z;   //!
   TBranch        *b_Particle_size;   //!
   TBranch        *b_GenVertex_;   //!
   TBranch        *b_GenVertex_fUniqueID;   //!
   TBranch        *b_GenVertex_fBits;   //!
   TBranch        *b_GenVertex_T;   //!
   TBranch        *b_GenVertex_X;   //!
   TBranch        *b_GenVertex_Y;   //!
   TBranch        *b_GenVertex_Z;   //!
   TBranch        *b_GenVertex_ErrorT;   //!
   TBranch        *b_GenVertex_ErrorX;   //!
   TBranch        *b_GenVertex_ErrorY;   //!
   TBranch        *b_GenVertex_ErrorZ;   //!
   TBranch        *b_GenVertex_Index;   //!
   TBranch        *b_GenVertex_NDF;   //!
   TBranch        *b_GenVertex_Sigma;   //!
   TBranch        *b_GenVertex_SumPT2;   //!
   TBranch        *b_GenVertex_GenSumPT2;   //!
   TBranch        *b_GenVertex_GenDeltaZ;   //!
   TBranch        *b_GenVertex_BTVSumPT2;   //!
   TBranch        *b_GenVertex_Constituents;   //!
   TBranch        *b_GenVertex_size;   //!
   TBranch        *b_Track_;   //!
   TBranch        *b_Track_fUniqueID;   //!
   TBranch        *b_Track_fBits;   //!
   TBranch        *b_Track_PID;   //!
   TBranch        *b_Track_Chi_k;   //!
   TBranch        *b_Track_Chi_pi;   //!
   TBranch        *b_Track_Counting_eff;   //!
   TBranch        *b_Track_L_DC;   //!
   TBranch        *b_Track_Truth_PID;   //!
   TBranch        *b_Track_Truth_P;   //!
   TBranch        *b_Track_Truth_CosTheta;   //!
   TBranch        *b_Track_Truth_Phi;   //!
   TBranch        *b_Track_Charge;   //!
   TBranch        *b_Track_P;   //!
   TBranch        *b_Track_PT;   //!
   TBranch        *b_Track_Eta;   //!
   TBranch        *b_Track_Phi;   //!
   TBranch        *b_Track_CtgTheta;   //!
   TBranch        *b_Track_CosTheta;   //!
   TBranch        *b_Track_C;   //!
   TBranch        *b_Track_Mass;   //!
   TBranch        *b_Track_Prob_K;   //!
   TBranch        *b_Track_Prob_Pi;   //!
   TBranch        *b_Track_Prob_P;   //!
   TBranch        *b_Track_EtaOuter;   //!
   TBranch        *b_Track_PhiOuter;   //!
   TBranch        *b_Track_T;   //!
   TBranch        *b_Track_X;   //!
   TBranch        *b_Track_Y;   //!
   TBranch        *b_Track_Z;   //!
   TBranch        *b_Track_TOuter;   //!
   TBranch        *b_Track_XOuter;   //!
   TBranch        *b_Track_YOuter;   //!
   TBranch        *b_Track_ZOuter;   //!
   TBranch        *b_Track_Xd;   //!
   TBranch        *b_Track_Yd;   //!
   TBranch        *b_Track_Zd;   //!
   TBranch        *b_Track_L;   //!
   TBranch        *b_Track_D0;   //!
   TBranch        *b_Track_DZ;   //!
   TBranch        *b_Track_Nclusters;   //!
   TBranch        *b_Track_Nclusters_err;   //!
   TBranch        *b_Track_dNdx;   //!
   TBranch        *b_Track_ErrorP;   //!
   TBranch        *b_Track_ErrorPT;   //!
   TBranch        *b_Track_ErrorPhi;   //!
   TBranch        *b_Track_ErrorCtgTheta;   //!
   TBranch        *b_Track_ErrorT;   //!
   TBranch        *b_Track_ErrorD0;   //!
   TBranch        *b_Track_ErrorDZ;   //!
   TBranch        *b_Track_ErrorC;   //!
   TBranch        *b_Track_ErrorD0Phi;   //!
   TBranch        *b_Track_ErrorD0C;   //!
   TBranch        *b_Track_ErrorD0DZ;   //!
   TBranch        *b_Track_ErrorD0CtgTheta;   //!
   TBranch        *b_Track_ErrorPhiC;   //!
   TBranch        *b_Track_ErrorPhiDZ;   //!
   TBranch        *b_Track_ErrorPhiCtgTheta;   //!
   TBranch        *b_Track_ErrorCDZ;   //!
   TBranch        *b_Track_ErrorCCtgTheta;   //!
   TBranch        *b_Track_ErrorDZCtgTheta;   //!
   TBranch        *b_Track_Particle;   //!
   TBranch        *b_Track_VertexIndex;   //!
   TBranch        *b_Track_size;   //!
   TBranch        *b_EFlowTrack_;   //!
   TBranch        *b_EFlowTrack_fUniqueID;   //!
   TBranch        *b_EFlowTrack_fBits;   //!
   TBranch        *b_EFlowTrack_PID;   //!
   TBranch        *b_EFlowTrack_Chi_k;   //!
   TBranch        *b_EFlowTrack_Chi_pi;   //!
   TBranch        *b_EFlowTrack_Counting_eff;   //!
   TBranch        *b_EFlowTrack_L_DC;   //!
   TBranch        *b_EFlowTrack_Truth_PID;   //!
   TBranch        *b_EFlowTrack_Truth_P;   //!
   TBranch        *b_EFlowTrack_Truth_CosTheta;   //!
   TBranch        *b_EFlowTrack_Truth_Phi;   //!
   TBranch        *b_EFlowTrack_Charge;   //!
   TBranch        *b_EFlowTrack_P;   //!
   TBranch        *b_EFlowTrack_PT;   //!
   TBranch        *b_EFlowTrack_Eta;   //!
   TBranch        *b_EFlowTrack_Phi;   //!
   TBranch        *b_EFlowTrack_CtgTheta;   //!
   TBranch        *b_EFlowTrack_CosTheta;   //!
   TBranch        *b_EFlowTrack_C;   //!
   TBranch        *b_EFlowTrack_Mass;   //!
   TBranch        *b_EFlowTrack_Prob_K;   //!
   TBranch        *b_EFlowTrack_Prob_Pi;   //!
   TBranch        *b_EFlowTrack_Prob_P;   //!
   TBranch        *b_EFlowTrack_EtaOuter;   //!
   TBranch        *b_EFlowTrack_PhiOuter;   //!
   TBranch        *b_EFlowTrack_T;   //!
   TBranch        *b_EFlowTrack_X;   //!
   TBranch        *b_EFlowTrack_Y;   //!
   TBranch        *b_EFlowTrack_Z;   //!
   TBranch        *b_EFlowTrack_TOuter;   //!
   TBranch        *b_EFlowTrack_XOuter;   //!
   TBranch        *b_EFlowTrack_YOuter;   //!
   TBranch        *b_EFlowTrack_ZOuter;   //!
   TBranch        *b_EFlowTrack_Xd;   //!
   TBranch        *b_EFlowTrack_Yd;   //!
   TBranch        *b_EFlowTrack_Zd;   //!
   TBranch        *b_EFlowTrack_L;   //!
   TBranch        *b_EFlowTrack_D0;   //!
   TBranch        *b_EFlowTrack_DZ;   //!
   TBranch        *b_EFlowTrack_Nclusters;   //!
   TBranch        *b_EFlowTrack_Nclusters_err;   //!
   TBranch        *b_EFlowTrack_dNdx;   //!
   TBranch        *b_EFlowTrack_ErrorP;   //!
   TBranch        *b_EFlowTrack_ErrorPT;   //!
   TBranch        *b_EFlowTrack_ErrorPhi;   //!
   TBranch        *b_EFlowTrack_ErrorCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorT;   //!
   TBranch        *b_EFlowTrack_ErrorD0;   //!
   TBranch        *b_EFlowTrack_ErrorDZ;   //!
   TBranch        *b_EFlowTrack_ErrorC;   //!
   TBranch        *b_EFlowTrack_ErrorD0Phi;   //!
   TBranch        *b_EFlowTrack_ErrorD0C;   //!
   TBranch        *b_EFlowTrack_ErrorD0DZ;   //!
   TBranch        *b_EFlowTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorPhiC;   //!
   TBranch        *b_EFlowTrack_ErrorPhiDZ;   //!
   TBranch        *b_EFlowTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorCDZ;   //!
   TBranch        *b_EFlowTrack_ErrorCCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_EFlowTrack_Particle;   //!
   TBranch        *b_EFlowTrack_VertexIndex;   //!
   TBranch        *b_EFlowTrack_size;   //!
   TBranch        *b_ParticleFlowCandidate_;   //!
   TBranch        *b_ParticleFlowCandidate_fUniqueID;   //!
   TBranch        *b_ParticleFlowCandidate_fBits;   //!
   TBranch        *b_ParticleFlowCandidate_PID;   //!
   TBranch        *b_ParticleFlowCandidate_Charge;   //!
   TBranch        *b_ParticleFlowCandidate_E;   //!
   TBranch        *b_ParticleFlowCandidate_P;   //!
   TBranch        *b_ParticleFlowCandidate_PT;   //!
   TBranch        *b_ParticleFlowCandidate_Eta;   //!
   TBranch        *b_ParticleFlowCandidate_Phi;   //!
   TBranch        *b_ParticleFlowCandidate_CtgTheta;   //!
   TBranch        *b_ParticleFlowCandidate_C;   //!
   TBranch        *b_ParticleFlowCandidate_Mass;   //!
   TBranch        *b_ParticleFlowCandidate_EtaOuter;   //!
   TBranch        *b_ParticleFlowCandidate_PhiOuter;   //!
   TBranch        *b_ParticleFlowCandidate_T;   //!
   TBranch        *b_ParticleFlowCandidate_X;   //!
   TBranch        *b_ParticleFlowCandidate_Y;   //!
   TBranch        *b_ParticleFlowCandidate_Z;   //!
   TBranch        *b_ParticleFlowCandidate_TOuter;   //!
   TBranch        *b_ParticleFlowCandidate_XOuter;   //!
   TBranch        *b_ParticleFlowCandidate_YOuter;   //!
   TBranch        *b_ParticleFlowCandidate_ZOuter;   //!
   TBranch        *b_ParticleFlowCandidate_Xd;   //!
   TBranch        *b_ParticleFlowCandidate_Yd;   //!
   TBranch        *b_ParticleFlowCandidate_Zd;   //!
   TBranch        *b_ParticleFlowCandidate_L;   //!
   TBranch        *b_ParticleFlowCandidate_D0;   //!
   TBranch        *b_ParticleFlowCandidate_DZ;   //!
   TBranch        *b_ParticleFlowCandidate_Nclusters;   //!
   TBranch        *b_ParticleFlowCandidate_Nclusters_err;   //!
   TBranch        *b_ParticleFlowCandidate_dNdx;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorP;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorPT;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorPhi;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorCtgTheta;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorT;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorD0;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorDZ;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorC;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorD0Phi;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorD0C;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorD0DZ;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorD0CtgTheta;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorPhiC;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorPhiDZ;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorPhiCtgTheta;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorCDZ;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorCCtgTheta;   //!
   TBranch        *b_ParticleFlowCandidate_ErrorDZCtgTheta;   //!
   TBranch        *b_ParticleFlowCandidate_VertexIndex;   //!
   TBranch        *b_ParticleFlowCandidate_NTimeHits;   //!
   TBranch        *b_ParticleFlowCandidate_Eem;   //!
   TBranch        *b_ParticleFlowCandidate_Ehad;   //!
   TBranch        *b_ParticleFlowCandidate_Etrk;   //!
   TBranch        *b_ParticleFlowCandidate_Edges;   //!
   TBranch        *b_ParticleFlowCandidate_Particles;   //!
   TBranch        *b_ParticleFlowCandidate_size;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_fUniqueID;   //!
   TBranch        *b_Electron_fBits;   //!
   TBranch        *b_Electron_PT;   //!
   TBranch        *b_Electron_Eta;   //!
   TBranch        *b_Electron_Phi;   //!
   TBranch        *b_Electron_T;   //!
   TBranch        *b_Electron_Charge;   //!
   TBranch        *b_Electron_EhadOverEem;   //!
   TBranch        *b_Electron_Particle;   //!
   TBranch        *b_Electron_IsolationVar;   //!
   TBranch        *b_Electron_IsolationVarRhoCorr;   //!
   TBranch        *b_Electron_SumPtCharged;   //!
   TBranch        *b_Electron_SumPtNeutral;   //!
   TBranch        *b_Electron_SumPtChargedPU;   //!
   TBranch        *b_Electron_SumPt;   //!
   TBranch        *b_Electron_D0;   //!
   TBranch        *b_Electron_DZ;   //!
   TBranch        *b_Electron_ErrorD0;   //!
   TBranch        *b_Electron_ErrorDZ;   //!
   TBranch        *b_Electron_size;   //!
   TBranch        *b_Photon_;   //!
   TBranch        *b_Photon_fUniqueID;   //!
   TBranch        *b_Photon_fBits;   //!
   TBranch        *b_Photon_PT;   //!
   TBranch        *b_Photon_Eta;   //!
   TBranch        *b_Photon_Phi;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_T;   //!
   TBranch        *b_Photon_EhadOverEem;   //!
   TBranch        *b_Photon_Particles;   //!
   TBranch        *b_Photon_IsolationVar;   //!
   TBranch        *b_Photon_IsolationVarRhoCorr;   //!
   TBranch        *b_Photon_SumPtCharged;   //!
   TBranch        *b_Photon_SumPtNeutral;   //!
   TBranch        *b_Photon_SumPtChargedPU;   //!
   TBranch        *b_Photon_SumPt;   //!
   TBranch        *b_Photon_CosTheta;   //!
   TBranch        *b_Photon_Truth_Phi;   //!
   TBranch        *b_Photon_Truth_E;   //!
   TBranch        *b_Photon_Truth_CosTheta;   //!
   TBranch        *b_Photon_Status;   //!
   TBranch        *b_Photon_size;   //!
   TBranch        *b_Muon_;   //!
   TBranch        *b_Muon_fUniqueID;   //!
   TBranch        *b_Muon_fBits;   //!
   TBranch        *b_Muon_PT;   //!
   TBranch        *b_Muon_Eta;   //!
   TBranch        *b_Muon_Phi;   //!
   TBranch        *b_Muon_T;   //!
   TBranch        *b_Muon_Charge;   //!
   TBranch        *b_Muon_Particle;   //!
   TBranch        *b_Muon_IsolationVar;   //!
   TBranch        *b_Muon_IsolationVarRhoCorr;   //!
   TBranch        *b_Muon_SumPtCharged;   //!
   TBranch        *b_Muon_SumPtNeutral;   //!
   TBranch        *b_Muon_SumPtChargedPU;   //!
   TBranch        *b_Muon_SumPt;   //!
   TBranch        *b_Muon_D0;   //!
   TBranch        *b_Muon_DZ;   //!
   TBranch        *b_Muon_ErrorD0;   //!
   TBranch        *b_Muon_ErrorDZ;   //!
   TBranch        *b_Muon_size;   //!

   b(TTree *tree=0);
   virtual ~b();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef b_cxx
b::b(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../rootfile/B.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../rootfile/B.root");
      }
      f->GetObject("Delphes",tree);

   }
   Init(tree);
}

b::~b()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t b::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t b::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void b::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
   fChain->SetBranchAddress("Event.Number", Event_Number, &b_Event_Number);
   fChain->SetBranchAddress("Event.ReadTime", Event_ReadTime, &b_Event_ReadTime);
   fChain->SetBranchAddress("Event.ProcTime", Event_ProcTime, &b_Event_ProcTime);
   fChain->SetBranchAddress("Event.ProcessID", Event_ProcessID, &b_Event_ProcessID);
   fChain->SetBranchAddress("Event.MPI", Event_MPI, &b_Event_MPI);
   fChain->SetBranchAddress("Event.Weight", Event_Weight, &b_Event_Weight);
   fChain->SetBranchAddress("Event.CrossSection", Event_CrossSection, &b_Event_CrossSection);
   fChain->SetBranchAddress("Event.CrossSectionError", Event_CrossSectionError, &b_Event_CrossSectionError);
   fChain->SetBranchAddress("Event.Scale", Event_Scale, &b_Event_Scale);
   fChain->SetBranchAddress("Event.AlphaQED", Event_AlphaQED, &b_Event_AlphaQED);
   fChain->SetBranchAddress("Event.AlphaQCD", Event_AlphaQCD, &b_Event_AlphaQCD);
   fChain->SetBranchAddress("Event.ID1", Event_ID1, &b_Event_ID1);
   fChain->SetBranchAddress("Event.ID2", Event_ID2, &b_Event_ID2);
   fChain->SetBranchAddress("Event.X1", Event_X1, &b_Event_X1);
   fChain->SetBranchAddress("Event.X2", Event_X2, &b_Event_X2);
   fChain->SetBranchAddress("Event.ScalePDF", Event_ScalePDF, &b_Event_ScalePDF);
   fChain->SetBranchAddress("Event.PDF1", Event_PDF1, &b_Event_PDF1);
   fChain->SetBranchAddress("Event.PDF2", Event_PDF2, &b_Event_PDF2);
   fChain->SetBranchAddress("Event_size", &Event_size, &b_Event_size);
   fChain->SetBranchAddress("Weight", &Weight_, &b_Weight_);
   fChain->SetBranchAddress("Weight.fUniqueID", Weight_fUniqueID, &b_Weight_fUniqueID);
   fChain->SetBranchAddress("Weight.fBits", Weight_fBits, &b_Weight_fBits);
   fChain->SetBranchAddress("Weight.Weight", Weight_Weight, &b_Weight_Weight);
   fChain->SetBranchAddress("Weight_size", &Weight_size, &b_Weight_size);
   fChain->SetBranchAddress("Particle", &Particle_, &b_Particle_);
   fChain->SetBranchAddress("Particle.fUniqueID", Particle_fUniqueID, &b_Particle_fUniqueID);
   fChain->SetBranchAddress("Particle.fBits", Particle_fBits, &b_Particle_fBits);
   fChain->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
   fChain->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
   fChain->SetBranchAddress("Particle.IsPU", Particle_IsPU, &b_Particle_IsPU);
   fChain->SetBranchAddress("Particle.M1", Particle_M1, &b_Particle_M1);
   fChain->SetBranchAddress("Particle.M2", Particle_M2, &b_Particle_M2);
   fChain->SetBranchAddress("Particle.D1", Particle_D1, &b_Particle_D1);
   fChain->SetBranchAddress("Particle.D2", Particle_D2, &b_Particle_D2);
   fChain->SetBranchAddress("Particle.Charge", Particle_Charge, &b_Particle_Charge);
   fChain->SetBranchAddress("Particle.Mass", Particle_Mass, &b_Particle_Mass);
   fChain->SetBranchAddress("Particle.E", Particle_E, &b_Particle_E);
   fChain->SetBranchAddress("Particle.Px", Particle_Px, &b_Particle_Px);
   fChain->SetBranchAddress("Particle.Py", Particle_Py, &b_Particle_Py);
   fChain->SetBranchAddress("Particle.Pz", Particle_Pz, &b_Particle_Pz);
   fChain->SetBranchAddress("Particle.P", Particle_P, &b_Particle_P);
   fChain->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
   fChain->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
   fChain->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
   fChain->SetBranchAddress("Particle.Rapidity", Particle_Rapidity, &b_Particle_Rapidity);
   fChain->SetBranchAddress("Particle.T", Particle_T, &b_Particle_T);
   fChain->SetBranchAddress("Particle.X", Particle_X, &b_Particle_X);
   fChain->SetBranchAddress("Particle.Y", Particle_Y, &b_Particle_Y);
   fChain->SetBranchAddress("Particle.Z", Particle_Z, &b_Particle_Z);
   fChain->SetBranchAddress("Particle_size", &Particle_size, &b_Particle_size);
   fChain->SetBranchAddress("GenVertex", &GenVertex_, &b_GenVertex_);
   fChain->SetBranchAddress("GenVertex.fUniqueID", GenVertex_fUniqueID, &b_GenVertex_fUniqueID);
   fChain->SetBranchAddress("GenVertex.fBits", GenVertex_fBits, &b_GenVertex_fBits);
   fChain->SetBranchAddress("GenVertex.T", GenVertex_T, &b_GenVertex_T);
   fChain->SetBranchAddress("GenVertex.X", GenVertex_X, &b_GenVertex_X);
   fChain->SetBranchAddress("GenVertex.Y", GenVertex_Y, &b_GenVertex_Y);
   fChain->SetBranchAddress("GenVertex.Z", GenVertex_Z, &b_GenVertex_Z);
   fChain->SetBranchAddress("GenVertex.ErrorT", GenVertex_ErrorT, &b_GenVertex_ErrorT);
   fChain->SetBranchAddress("GenVertex.ErrorX", GenVertex_ErrorX, &b_GenVertex_ErrorX);
   fChain->SetBranchAddress("GenVertex.ErrorY", GenVertex_ErrorY, &b_GenVertex_ErrorY);
   fChain->SetBranchAddress("GenVertex.ErrorZ", GenVertex_ErrorZ, &b_GenVertex_ErrorZ);
   fChain->SetBranchAddress("GenVertex.Index", GenVertex_Index, &b_GenVertex_Index);
   fChain->SetBranchAddress("GenVertex.NDF", GenVertex_NDF, &b_GenVertex_NDF);
   fChain->SetBranchAddress("GenVertex.Sigma", GenVertex_Sigma, &b_GenVertex_Sigma);
   fChain->SetBranchAddress("GenVertex.SumPT2", GenVertex_SumPT2, &b_GenVertex_SumPT2);
   fChain->SetBranchAddress("GenVertex.GenSumPT2", GenVertex_GenSumPT2, &b_GenVertex_GenSumPT2);
   fChain->SetBranchAddress("GenVertex.GenDeltaZ", GenVertex_GenDeltaZ, &b_GenVertex_GenDeltaZ);
   fChain->SetBranchAddress("GenVertex.BTVSumPT2", GenVertex_BTVSumPT2, &b_GenVertex_BTVSumPT2);
   fChain->SetBranchAddress("GenVertex.Constituents", GenVertex_Constituents, &b_GenVertex_Constituents);
   fChain->SetBranchAddress("GenVertex_size", &GenVertex_size, &b_GenVertex_size);
   fChain->SetBranchAddress("Track", &Track_, &b_Track_);
   fChain->SetBranchAddress("Track.fUniqueID", Track_fUniqueID, &b_Track_fUniqueID);
   fChain->SetBranchAddress("Track.fBits", Track_fBits, &b_Track_fBits);
   fChain->SetBranchAddress("Track.PID", Track_PID, &b_Track_PID);
   fChain->SetBranchAddress("Track.Chi_k", Track_Chi_k, &b_Track_Chi_k);
   fChain->SetBranchAddress("Track.Chi_pi", Track_Chi_pi, &b_Track_Chi_pi);
   fChain->SetBranchAddress("Track.Counting_eff", Track_Counting_eff, &b_Track_Counting_eff);
   fChain->SetBranchAddress("Track.L_DC", Track_L_DC, &b_Track_L_DC);
   fChain->SetBranchAddress("Track.Truth_PID", Track_Truth_PID, &b_Track_Truth_PID);
   fChain->SetBranchAddress("Track.Truth_P", Track_Truth_P, &b_Track_Truth_P);
   fChain->SetBranchAddress("Track.Truth_CosTheta", Track_Truth_CosTheta, &b_Track_Truth_CosTheta);
   fChain->SetBranchAddress("Track.Truth_Phi", Track_Truth_Phi, &b_Track_Truth_Phi);
   fChain->SetBranchAddress("Track.Charge", Track_Charge, &b_Track_Charge);
   fChain->SetBranchAddress("Track.P", Track_P, &b_Track_P);
   fChain->SetBranchAddress("Track.PT", Track_PT, &b_Track_PT);
   fChain->SetBranchAddress("Track.Eta", Track_Eta, &b_Track_Eta);
   fChain->SetBranchAddress("Track.Phi", Track_Phi, &b_Track_Phi);
   fChain->SetBranchAddress("Track.CtgTheta", Track_CtgTheta, &b_Track_CtgTheta);
   fChain->SetBranchAddress("Track.CosTheta", Track_CosTheta, &b_Track_CosTheta);
   fChain->SetBranchAddress("Track.C", Track_C, &b_Track_C);
   fChain->SetBranchAddress("Track.Mass", Track_Mass, &b_Track_Mass);
   fChain->SetBranchAddress("Track.Prob_K", Track_Prob_K, &b_Track_Prob_K);
   fChain->SetBranchAddress("Track.Prob_Pi", Track_Prob_Pi, &b_Track_Prob_Pi);
   fChain->SetBranchAddress("Track.Prob_P", Track_Prob_P, &b_Track_Prob_P);
   fChain->SetBranchAddress("Track.EtaOuter", Track_EtaOuter, &b_Track_EtaOuter);
   fChain->SetBranchAddress("Track.PhiOuter", Track_PhiOuter, &b_Track_PhiOuter);
   fChain->SetBranchAddress("Track.T", Track_T, &b_Track_T);
   fChain->SetBranchAddress("Track.X", Track_X, &b_Track_X);
   fChain->SetBranchAddress("Track.Y", Track_Y, &b_Track_Y);
   fChain->SetBranchAddress("Track.Z", Track_Z, &b_Track_Z);
   fChain->SetBranchAddress("Track.TOuter", Track_TOuter, &b_Track_TOuter);
   fChain->SetBranchAddress("Track.XOuter", Track_XOuter, &b_Track_XOuter);
   fChain->SetBranchAddress("Track.YOuter", Track_YOuter, &b_Track_YOuter);
   fChain->SetBranchAddress("Track.ZOuter", Track_ZOuter, &b_Track_ZOuter);
   fChain->SetBranchAddress("Track.Xd", Track_Xd, &b_Track_Xd);
   fChain->SetBranchAddress("Track.Yd", Track_Yd, &b_Track_Yd);
   fChain->SetBranchAddress("Track.Zd", Track_Zd, &b_Track_Zd);
   fChain->SetBranchAddress("Track.L", Track_L, &b_Track_L);
   fChain->SetBranchAddress("Track.D0", Track_D0, &b_Track_D0);
   fChain->SetBranchAddress("Track.DZ", Track_DZ, &b_Track_DZ);
   fChain->SetBranchAddress("Track.Nclusters", Track_Nclusters, &b_Track_Nclusters);
   fChain->SetBranchAddress("Track.Nclusters_err", Track_Nclusters_err, &b_Track_Nclusters_err);
   fChain->SetBranchAddress("Track.dNdx", Track_dNdx, &b_Track_dNdx);
   fChain->SetBranchAddress("Track.ErrorP", Track_ErrorP, &b_Track_ErrorP);
   fChain->SetBranchAddress("Track.ErrorPT", Track_ErrorPT, &b_Track_ErrorPT);
   fChain->SetBranchAddress("Track.ErrorPhi", Track_ErrorPhi, &b_Track_ErrorPhi);
   fChain->SetBranchAddress("Track.ErrorCtgTheta", Track_ErrorCtgTheta, &b_Track_ErrorCtgTheta);
   fChain->SetBranchAddress("Track.ErrorT", Track_ErrorT, &b_Track_ErrorT);
   fChain->SetBranchAddress("Track.ErrorD0", Track_ErrorD0, &b_Track_ErrorD0);
   fChain->SetBranchAddress("Track.ErrorDZ", Track_ErrorDZ, &b_Track_ErrorDZ);
   fChain->SetBranchAddress("Track.ErrorC", Track_ErrorC, &b_Track_ErrorC);
   fChain->SetBranchAddress("Track.ErrorD0Phi", Track_ErrorD0Phi, &b_Track_ErrorD0Phi);
   fChain->SetBranchAddress("Track.ErrorD0C", Track_ErrorD0C, &b_Track_ErrorD0C);
   fChain->SetBranchAddress("Track.ErrorD0DZ", Track_ErrorD0DZ, &b_Track_ErrorD0DZ);
   fChain->SetBranchAddress("Track.ErrorD0CtgTheta", Track_ErrorD0CtgTheta, &b_Track_ErrorD0CtgTheta);
   fChain->SetBranchAddress("Track.ErrorPhiC", Track_ErrorPhiC, &b_Track_ErrorPhiC);
   fChain->SetBranchAddress("Track.ErrorPhiDZ", Track_ErrorPhiDZ, &b_Track_ErrorPhiDZ);
   fChain->SetBranchAddress("Track.ErrorPhiCtgTheta", Track_ErrorPhiCtgTheta, &b_Track_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("Track.ErrorCDZ", Track_ErrorCDZ, &b_Track_ErrorCDZ);
   fChain->SetBranchAddress("Track.ErrorCCtgTheta", Track_ErrorCCtgTheta, &b_Track_ErrorCCtgTheta);
   fChain->SetBranchAddress("Track.ErrorDZCtgTheta", Track_ErrorDZCtgTheta, &b_Track_ErrorDZCtgTheta);
   fChain->SetBranchAddress("Track.Particle", Track_Particle, &b_Track_Particle);
   fChain->SetBranchAddress("Track.VertexIndex", Track_VertexIndex, &b_Track_VertexIndex);
   fChain->SetBranchAddress("Track_size", &Track_size, &b_Track_size);
   fChain->SetBranchAddress("EFlowTrack", &EFlowTrack_, &b_EFlowTrack_);
   fChain->SetBranchAddress("EFlowTrack.fUniqueID", EFlowTrack_fUniqueID, &b_EFlowTrack_fUniqueID);
   fChain->SetBranchAddress("EFlowTrack.fBits", EFlowTrack_fBits, &b_EFlowTrack_fBits);
   fChain->SetBranchAddress("EFlowTrack.PID", EFlowTrack_PID, &b_EFlowTrack_PID);
   fChain->SetBranchAddress("EFlowTrack.Chi_k", EFlowTrack_Chi_k, &b_EFlowTrack_Chi_k);
   fChain->SetBranchAddress("EFlowTrack.Chi_pi", EFlowTrack_Chi_pi, &b_EFlowTrack_Chi_pi);
   fChain->SetBranchAddress("EFlowTrack.Counting_eff", EFlowTrack_Counting_eff, &b_EFlowTrack_Counting_eff);
   fChain->SetBranchAddress("EFlowTrack.L_DC", EFlowTrack_L_DC, &b_EFlowTrack_L_DC);
   fChain->SetBranchAddress("EFlowTrack.Truth_PID", EFlowTrack_Truth_PID, &b_EFlowTrack_Truth_PID);
   fChain->SetBranchAddress("EFlowTrack.Truth_P", EFlowTrack_Truth_P, &b_EFlowTrack_Truth_P);
   fChain->SetBranchAddress("EFlowTrack.Truth_CosTheta", EFlowTrack_Truth_CosTheta, &b_EFlowTrack_Truth_CosTheta);
   fChain->SetBranchAddress("EFlowTrack.Truth_Phi", EFlowTrack_Truth_Phi, &b_EFlowTrack_Truth_Phi);
   fChain->SetBranchAddress("EFlowTrack.Charge", EFlowTrack_Charge, &b_EFlowTrack_Charge);
   fChain->SetBranchAddress("EFlowTrack.P", EFlowTrack_P, &b_EFlowTrack_P);
   fChain->SetBranchAddress("EFlowTrack.PT", EFlowTrack_PT, &b_EFlowTrack_PT);
   fChain->SetBranchAddress("EFlowTrack.Eta", EFlowTrack_Eta, &b_EFlowTrack_Eta);
   fChain->SetBranchAddress("EFlowTrack.Phi", EFlowTrack_Phi, &b_EFlowTrack_Phi);
   fChain->SetBranchAddress("EFlowTrack.CtgTheta", EFlowTrack_CtgTheta, &b_EFlowTrack_CtgTheta);
   fChain->SetBranchAddress("EFlowTrack.CosTheta", EFlowTrack_CosTheta, &b_EFlowTrack_CosTheta);
   fChain->SetBranchAddress("EFlowTrack.C", EFlowTrack_C, &b_EFlowTrack_C);
   fChain->SetBranchAddress("EFlowTrack.Mass", EFlowTrack_Mass, &b_EFlowTrack_Mass);
   fChain->SetBranchAddress("EFlowTrack.Prob_K", EFlowTrack_Prob_K, &b_EFlowTrack_Prob_K);
   fChain->SetBranchAddress("EFlowTrack.Prob_Pi", EFlowTrack_Prob_Pi, &b_EFlowTrack_Prob_Pi);
   fChain->SetBranchAddress("EFlowTrack.Prob_P", EFlowTrack_Prob_P, &b_EFlowTrack_Prob_P);
   fChain->SetBranchAddress("EFlowTrack.EtaOuter", EFlowTrack_EtaOuter, &b_EFlowTrack_EtaOuter);
   fChain->SetBranchAddress("EFlowTrack.PhiOuter", EFlowTrack_PhiOuter, &b_EFlowTrack_PhiOuter);
   fChain->SetBranchAddress("EFlowTrack.T", EFlowTrack_T, &b_EFlowTrack_T);
   fChain->SetBranchAddress("EFlowTrack.X", EFlowTrack_X, &b_EFlowTrack_X);
   fChain->SetBranchAddress("EFlowTrack.Y", EFlowTrack_Y, &b_EFlowTrack_Y);
   fChain->SetBranchAddress("EFlowTrack.Z", EFlowTrack_Z, &b_EFlowTrack_Z);
   fChain->SetBranchAddress("EFlowTrack.TOuter", EFlowTrack_TOuter, &b_EFlowTrack_TOuter);
   fChain->SetBranchAddress("EFlowTrack.XOuter", EFlowTrack_XOuter, &b_EFlowTrack_XOuter);
   fChain->SetBranchAddress("EFlowTrack.YOuter", EFlowTrack_YOuter, &b_EFlowTrack_YOuter);
   fChain->SetBranchAddress("EFlowTrack.ZOuter", EFlowTrack_ZOuter, &b_EFlowTrack_ZOuter);
   fChain->SetBranchAddress("EFlowTrack.Xd", EFlowTrack_Xd, &b_EFlowTrack_Xd);
   fChain->SetBranchAddress("EFlowTrack.Yd", EFlowTrack_Yd, &b_EFlowTrack_Yd);
   fChain->SetBranchAddress("EFlowTrack.Zd", EFlowTrack_Zd, &b_EFlowTrack_Zd);
   fChain->SetBranchAddress("EFlowTrack.L", EFlowTrack_L, &b_EFlowTrack_L);
   fChain->SetBranchAddress("EFlowTrack.D0", EFlowTrack_D0, &b_EFlowTrack_D0);
   fChain->SetBranchAddress("EFlowTrack.DZ", EFlowTrack_DZ, &b_EFlowTrack_DZ);
   fChain->SetBranchAddress("EFlowTrack.Nclusters", EFlowTrack_Nclusters, &b_EFlowTrack_Nclusters);
   fChain->SetBranchAddress("EFlowTrack.Nclusters_err", EFlowTrack_Nclusters_err, &b_EFlowTrack_Nclusters_err);
   fChain->SetBranchAddress("EFlowTrack.dNdx", EFlowTrack_dNdx, &b_EFlowTrack_dNdx);
   fChain->SetBranchAddress("EFlowTrack.ErrorP", EFlowTrack_ErrorP, &b_EFlowTrack_ErrorP);
   fChain->SetBranchAddress("EFlowTrack.ErrorPT", EFlowTrack_ErrorPT, &b_EFlowTrack_ErrorPT);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhi", EFlowTrack_ErrorPhi, &b_EFlowTrack_ErrorPhi);
   fChain->SetBranchAddress("EFlowTrack.ErrorCtgTheta", EFlowTrack_ErrorCtgTheta, &b_EFlowTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorT", EFlowTrack_ErrorT, &b_EFlowTrack_ErrorT);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0", EFlowTrack_ErrorD0, &b_EFlowTrack_ErrorD0);
   fChain->SetBranchAddress("EFlowTrack.ErrorDZ", EFlowTrack_ErrorDZ, &b_EFlowTrack_ErrorDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorC", EFlowTrack_ErrorC, &b_EFlowTrack_ErrorC);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0Phi", EFlowTrack_ErrorD0Phi, &b_EFlowTrack_ErrorD0Phi);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0C", EFlowTrack_ErrorD0C, &b_EFlowTrack_ErrorD0C);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0DZ", EFlowTrack_ErrorD0DZ, &b_EFlowTrack_ErrorD0DZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0CtgTheta", EFlowTrack_ErrorD0CtgTheta, &b_EFlowTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiC", EFlowTrack_ErrorPhiC, &b_EFlowTrack_ErrorPhiC);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiDZ", EFlowTrack_ErrorPhiDZ, &b_EFlowTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiCtgTheta", EFlowTrack_ErrorPhiCtgTheta, &b_EFlowTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorCDZ", EFlowTrack_ErrorCDZ, &b_EFlowTrack_ErrorCDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorCCtgTheta", EFlowTrack_ErrorCCtgTheta, &b_EFlowTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorDZCtgTheta", EFlowTrack_ErrorDZCtgTheta, &b_EFlowTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.Particle", EFlowTrack_Particle, &b_EFlowTrack_Particle);
   fChain->SetBranchAddress("EFlowTrack.VertexIndex", EFlowTrack_VertexIndex, &b_EFlowTrack_VertexIndex);
   fChain->SetBranchAddress("EFlowTrack_size", &EFlowTrack_size, &b_EFlowTrack_size);
   fChain->SetBranchAddress("ParticleFlowCandidate", &ParticleFlowCandidate_, &b_ParticleFlowCandidate_);
   fChain->SetBranchAddress("ParticleFlowCandidate.fUniqueID", ParticleFlowCandidate_fUniqueID, &b_ParticleFlowCandidate_fUniqueID);
   fChain->SetBranchAddress("ParticleFlowCandidate.fBits", ParticleFlowCandidate_fBits, &b_ParticleFlowCandidate_fBits);
   fChain->SetBranchAddress("ParticleFlowCandidate.PID", ParticleFlowCandidate_PID, &b_ParticleFlowCandidate_PID);
   fChain->SetBranchAddress("ParticleFlowCandidate.Charge", ParticleFlowCandidate_Charge, &b_ParticleFlowCandidate_Charge);
   fChain->SetBranchAddress("ParticleFlowCandidate.E", ParticleFlowCandidate_E, &b_ParticleFlowCandidate_E);
   fChain->SetBranchAddress("ParticleFlowCandidate.P", ParticleFlowCandidate_P, &b_ParticleFlowCandidate_P);
   fChain->SetBranchAddress("ParticleFlowCandidate.PT", ParticleFlowCandidate_PT, &b_ParticleFlowCandidate_PT);
   fChain->SetBranchAddress("ParticleFlowCandidate.Eta", ParticleFlowCandidate_Eta, &b_ParticleFlowCandidate_Eta);
   fChain->SetBranchAddress("ParticleFlowCandidate.Phi", ParticleFlowCandidate_Phi, &b_ParticleFlowCandidate_Phi);
   fChain->SetBranchAddress("ParticleFlowCandidate.CtgTheta", ParticleFlowCandidate_CtgTheta, &b_ParticleFlowCandidate_CtgTheta);
   fChain->SetBranchAddress("ParticleFlowCandidate.C", ParticleFlowCandidate_C, &b_ParticleFlowCandidate_C);
   fChain->SetBranchAddress("ParticleFlowCandidate.Mass", ParticleFlowCandidate_Mass, &b_ParticleFlowCandidate_Mass);
   fChain->SetBranchAddress("ParticleFlowCandidate.EtaOuter", ParticleFlowCandidate_EtaOuter, &b_ParticleFlowCandidate_EtaOuter);
   fChain->SetBranchAddress("ParticleFlowCandidate.PhiOuter", ParticleFlowCandidate_PhiOuter, &b_ParticleFlowCandidate_PhiOuter);
   fChain->SetBranchAddress("ParticleFlowCandidate.T", ParticleFlowCandidate_T, &b_ParticleFlowCandidate_T);
   fChain->SetBranchAddress("ParticleFlowCandidate.X", ParticleFlowCandidate_X, &b_ParticleFlowCandidate_X);
   fChain->SetBranchAddress("ParticleFlowCandidate.Y", ParticleFlowCandidate_Y, &b_ParticleFlowCandidate_Y);
   fChain->SetBranchAddress("ParticleFlowCandidate.Z", ParticleFlowCandidate_Z, &b_ParticleFlowCandidate_Z);
   fChain->SetBranchAddress("ParticleFlowCandidate.TOuter", ParticleFlowCandidate_TOuter, &b_ParticleFlowCandidate_TOuter);
   fChain->SetBranchAddress("ParticleFlowCandidate.XOuter", ParticleFlowCandidate_XOuter, &b_ParticleFlowCandidate_XOuter);
   fChain->SetBranchAddress("ParticleFlowCandidate.YOuter", ParticleFlowCandidate_YOuter, &b_ParticleFlowCandidate_YOuter);
   fChain->SetBranchAddress("ParticleFlowCandidate.ZOuter", ParticleFlowCandidate_ZOuter, &b_ParticleFlowCandidate_ZOuter);
   fChain->SetBranchAddress("ParticleFlowCandidate.Xd", ParticleFlowCandidate_Xd, &b_ParticleFlowCandidate_Xd);
   fChain->SetBranchAddress("ParticleFlowCandidate.Yd", ParticleFlowCandidate_Yd, &b_ParticleFlowCandidate_Yd);
   fChain->SetBranchAddress("ParticleFlowCandidate.Zd", ParticleFlowCandidate_Zd, &b_ParticleFlowCandidate_Zd);
   fChain->SetBranchAddress("ParticleFlowCandidate.L", ParticleFlowCandidate_L, &b_ParticleFlowCandidate_L);
   fChain->SetBranchAddress("ParticleFlowCandidate.D0", ParticleFlowCandidate_D0, &b_ParticleFlowCandidate_D0);
   fChain->SetBranchAddress("ParticleFlowCandidate.DZ", ParticleFlowCandidate_DZ, &b_ParticleFlowCandidate_DZ);
   fChain->SetBranchAddress("ParticleFlowCandidate.Nclusters", ParticleFlowCandidate_Nclusters, &b_ParticleFlowCandidate_Nclusters);
   fChain->SetBranchAddress("ParticleFlowCandidate.Nclusters_err", ParticleFlowCandidate_Nclusters_err, &b_ParticleFlowCandidate_Nclusters_err);
   fChain->SetBranchAddress("ParticleFlowCandidate.dNdx", ParticleFlowCandidate_dNdx, &b_ParticleFlowCandidate_dNdx);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorP", ParticleFlowCandidate_ErrorP, &b_ParticleFlowCandidate_ErrorP);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorPT", ParticleFlowCandidate_ErrorPT, &b_ParticleFlowCandidate_ErrorPT);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorPhi", ParticleFlowCandidate_ErrorPhi, &b_ParticleFlowCandidate_ErrorPhi);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorCtgTheta", ParticleFlowCandidate_ErrorCtgTheta, &b_ParticleFlowCandidate_ErrorCtgTheta);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorT", ParticleFlowCandidate_ErrorT, &b_ParticleFlowCandidate_ErrorT);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorD0", ParticleFlowCandidate_ErrorD0, &b_ParticleFlowCandidate_ErrorD0);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorDZ", ParticleFlowCandidate_ErrorDZ, &b_ParticleFlowCandidate_ErrorDZ);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorC", ParticleFlowCandidate_ErrorC, &b_ParticleFlowCandidate_ErrorC);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorD0Phi", ParticleFlowCandidate_ErrorD0Phi, &b_ParticleFlowCandidate_ErrorD0Phi);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorD0C", ParticleFlowCandidate_ErrorD0C, &b_ParticleFlowCandidate_ErrorD0C);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorD0DZ", ParticleFlowCandidate_ErrorD0DZ, &b_ParticleFlowCandidate_ErrorD0DZ);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorD0CtgTheta", ParticleFlowCandidate_ErrorD0CtgTheta, &b_ParticleFlowCandidate_ErrorD0CtgTheta);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorPhiC", ParticleFlowCandidate_ErrorPhiC, &b_ParticleFlowCandidate_ErrorPhiC);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorPhiDZ", ParticleFlowCandidate_ErrorPhiDZ, &b_ParticleFlowCandidate_ErrorPhiDZ);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorPhiCtgTheta", ParticleFlowCandidate_ErrorPhiCtgTheta, &b_ParticleFlowCandidate_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorCDZ", ParticleFlowCandidate_ErrorCDZ, &b_ParticleFlowCandidate_ErrorCDZ);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorCCtgTheta", ParticleFlowCandidate_ErrorCCtgTheta, &b_ParticleFlowCandidate_ErrorCCtgTheta);
   fChain->SetBranchAddress("ParticleFlowCandidate.ErrorDZCtgTheta", ParticleFlowCandidate_ErrorDZCtgTheta, &b_ParticleFlowCandidate_ErrorDZCtgTheta);
   fChain->SetBranchAddress("ParticleFlowCandidate.VertexIndex", ParticleFlowCandidate_VertexIndex, &b_ParticleFlowCandidate_VertexIndex);
   fChain->SetBranchAddress("ParticleFlowCandidate.NTimeHits", ParticleFlowCandidate_NTimeHits, &b_ParticleFlowCandidate_NTimeHits);
   fChain->SetBranchAddress("ParticleFlowCandidate.Eem", ParticleFlowCandidate_Eem, &b_ParticleFlowCandidate_Eem);
   fChain->SetBranchAddress("ParticleFlowCandidate.Ehad", ParticleFlowCandidate_Ehad, &b_ParticleFlowCandidate_Ehad);
   fChain->SetBranchAddress("ParticleFlowCandidate.Etrk", ParticleFlowCandidate_Etrk, &b_ParticleFlowCandidate_Etrk);
   fChain->SetBranchAddress("ParticleFlowCandidate.Edges[4]", ParticleFlowCandidate_Edges, &b_ParticleFlowCandidate_Edges);
   fChain->SetBranchAddress("ParticleFlowCandidate.Particles", ParticleFlowCandidate_Particles, &b_ParticleFlowCandidate_Particles);
   fChain->SetBranchAddress("ParticleFlowCandidate_size", &ParticleFlowCandidate_size, &b_ParticleFlowCandidate_size);
   fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
   fChain->SetBranchAddress("Electron.fUniqueID", Electron_fUniqueID, &b_Electron_fUniqueID);
   fChain->SetBranchAddress("Electron.fBits", Electron_fBits, &b_Electron_fBits);
   fChain->SetBranchAddress("Electron.PT", Electron_PT, &b_Electron_PT);
   fChain->SetBranchAddress("Electron.Eta", Electron_Eta, &b_Electron_Eta);
   fChain->SetBranchAddress("Electron.Phi", Electron_Phi, &b_Electron_Phi);
   fChain->SetBranchAddress("Electron.T", Electron_T, &b_Electron_T);
   fChain->SetBranchAddress("Electron.Charge", Electron_Charge, &b_Electron_Charge);
   fChain->SetBranchAddress("Electron.EhadOverEem", Electron_EhadOverEem, &b_Electron_EhadOverEem);
   fChain->SetBranchAddress("Electron.Particle", Electron_Particle, &b_Electron_Particle);
   fChain->SetBranchAddress("Electron.IsolationVar", Electron_IsolationVar, &b_Electron_IsolationVar);
   fChain->SetBranchAddress("Electron.IsolationVarRhoCorr", Electron_IsolationVarRhoCorr, &b_Electron_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Electron.SumPtCharged", Electron_SumPtCharged, &b_Electron_SumPtCharged);
   fChain->SetBranchAddress("Electron.SumPtNeutral", Electron_SumPtNeutral, &b_Electron_SumPtNeutral);
   fChain->SetBranchAddress("Electron.SumPtChargedPU", Electron_SumPtChargedPU, &b_Electron_SumPtChargedPU);
   fChain->SetBranchAddress("Electron.SumPt", Electron_SumPt, &b_Electron_SumPt);
   fChain->SetBranchAddress("Electron.D0", Electron_D0, &b_Electron_D0);
   fChain->SetBranchAddress("Electron.DZ", Electron_DZ, &b_Electron_DZ);
   fChain->SetBranchAddress("Electron.ErrorD0", Electron_ErrorD0, &b_Electron_ErrorD0);
   fChain->SetBranchAddress("Electron.ErrorDZ", Electron_ErrorDZ, &b_Electron_ErrorDZ);
   fChain->SetBranchAddress("Electron_size", &Electron_size, &b_Electron_size);
   fChain->SetBranchAddress("Photon", &Photon_, &b_Photon_);
   fChain->SetBranchAddress("Photon.fUniqueID", Photon_fUniqueID, &b_Photon_fUniqueID);
   fChain->SetBranchAddress("Photon.fBits", Photon_fBits, &b_Photon_fBits);
   fChain->SetBranchAddress("Photon.PT", Photon_PT, &b_Photon_PT);
   fChain->SetBranchAddress("Photon.Eta", Photon_Eta, &b_Photon_Eta);
   fChain->SetBranchAddress("Photon.Phi", Photon_Phi, &b_Photon_Phi);
   fChain->SetBranchAddress("Photon.E", Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon.T", Photon_T, &b_Photon_T);
   fChain->SetBranchAddress("Photon.EhadOverEem", Photon_EhadOverEem, &b_Photon_EhadOverEem);
   fChain->SetBranchAddress("Photon.Particles", Photon_Particles, &b_Photon_Particles);
   fChain->SetBranchAddress("Photon.IsolationVar", Photon_IsolationVar, &b_Photon_IsolationVar);
   fChain->SetBranchAddress("Photon.IsolationVarRhoCorr", Photon_IsolationVarRhoCorr, &b_Photon_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Photon.SumPtCharged", Photon_SumPtCharged, &b_Photon_SumPtCharged);
   fChain->SetBranchAddress("Photon.SumPtNeutral", Photon_SumPtNeutral, &b_Photon_SumPtNeutral);
   fChain->SetBranchAddress("Photon.SumPtChargedPU", Photon_SumPtChargedPU, &b_Photon_SumPtChargedPU);
   fChain->SetBranchAddress("Photon.SumPt", Photon_SumPt, &b_Photon_SumPt);
   fChain->SetBranchAddress("Photon.CosTheta", Photon_CosTheta, &b_Photon_CosTheta);
   fChain->SetBranchAddress("Photon.Truth_Phi", Photon_Truth_Phi, &b_Photon_Truth_Phi);
   fChain->SetBranchAddress("Photon.Truth_E", Photon_Truth_E, &b_Photon_Truth_E);
   fChain->SetBranchAddress("Photon.Truth_CosTheta", Photon_Truth_CosTheta, &b_Photon_Truth_CosTheta);
   fChain->SetBranchAddress("Photon.Status", Photon_Status, &b_Photon_Status);
   fChain->SetBranchAddress("Photon_size", &Photon_size, &b_Photon_size);
   fChain->SetBranchAddress("Muon", &Muon_, &b_Muon_);
   fChain->SetBranchAddress("Muon.fUniqueID", Muon_fUniqueID, &b_Muon_fUniqueID);
   fChain->SetBranchAddress("Muon.fBits", Muon_fBits, &b_Muon_fBits);
   fChain->SetBranchAddress("Muon.PT", Muon_PT, &b_Muon_PT);
   fChain->SetBranchAddress("Muon.Eta", Muon_Eta, &b_Muon_Eta);
   fChain->SetBranchAddress("Muon.Phi", Muon_Phi, &b_Muon_Phi);
   fChain->SetBranchAddress("Muon.T", Muon_T, &b_Muon_T);
   fChain->SetBranchAddress("Muon.Charge", Muon_Charge, &b_Muon_Charge);
   fChain->SetBranchAddress("Muon.Particle", Muon_Particle, &b_Muon_Particle);
   fChain->SetBranchAddress("Muon.IsolationVar", Muon_IsolationVar, &b_Muon_IsolationVar);
   fChain->SetBranchAddress("Muon.IsolationVarRhoCorr", Muon_IsolationVarRhoCorr, &b_Muon_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Muon.SumPtCharged", Muon_SumPtCharged, &b_Muon_SumPtCharged);
   fChain->SetBranchAddress("Muon.SumPtNeutral", Muon_SumPtNeutral, &b_Muon_SumPtNeutral);
   fChain->SetBranchAddress("Muon.SumPtChargedPU", Muon_SumPtChargedPU, &b_Muon_SumPtChargedPU);
   fChain->SetBranchAddress("Muon.SumPt", Muon_SumPt, &b_Muon_SumPt);
   fChain->SetBranchAddress("Muon.D0", Muon_D0, &b_Muon_D0);
   fChain->SetBranchAddress("Muon.DZ", Muon_DZ, &b_Muon_DZ);
   fChain->SetBranchAddress("Muon.ErrorD0", Muon_ErrorD0, &b_Muon_ErrorD0);
   fChain->SetBranchAddress("Muon.ErrorDZ", Muon_ErrorDZ, &b_Muon_ErrorDZ);
   fChain->SetBranchAddress("Muon_size", &Muon_size, &b_Muon_size);
   Notify();
}

Bool_t b::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void b::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t b::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef b_cxx
