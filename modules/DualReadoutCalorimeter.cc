/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class DualReadoutCalorimeter
 *

  // ============  TODO  =========================================================
  // This implementation of dual calorimetry relies on several approximations:
  // - If hadronic energy is found in the tower the energy resolution then the full tower enrgy is smeared according to hadronic resolution (pessimistic for (e,n) or (pi+,gamma))
  // - While e+ vs pi+ (or gamma vs n) separation is in principle possible for single particles (using C/S, PMT timing, lateral shower profile) it is not obvious it can be done overlapping particles.
  //   Now we assume that regarless of the number of particle hits per tower we can always distinguish e+ vs pi+, which is probably not true in the case (e+,n) vs (pi+,gamma) without longitudinal segmentation.

 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "modules/DualReadoutCalorimeter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

DualReadoutCalorimeter::DualReadoutCalorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0)
{

  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;

  fECalTowerTrackArray = new TObjArray;
  fItECalTowerTrackArray = fECalTowerTrackArray->MakeIterator();

  fHCalTowerTrackArray = new TObjArray;
  fItHCalTowerTrackArray = fHCalTowerTrackArray->MakeIterator();

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();

}

//------------------------------------------------------------------------------

DualReadoutCalorimeter::~DualReadoutCalorimeter()
{

  if(fECalResolutionFormula) delete fECalResolutionFormula;
  if(fHCalResolutionFormula) delete fHCalResolutionFormula;

  if(fECalTowerTrackArray) delete fECalTowerTrackArray;
  if(fItECalTowerTrackArray) delete fItECalTowerTrackArray;

  if(fHCalTowerTrackArray) delete fHCalTowerTrackArray;
  if(fItHCalTowerTrackArray) delete fItHCalTowerTrackArray;

  if(fTowerTrackArray) delete fTowerTrackArray;
  if(fItTowerTrackArray) delete fItTowerTrackArray;

}

//------------------------------------------------------------------------------

void DualReadoutCalorimeter::Init()
{
  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  Double_t ecalFraction, hcalFraction;
  TBinMap::iterator itEtaBin;
  set< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size/2; ++i)
  {
    paramEtaBins = param[i*2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i*2 + 1];
    sizePhiBins = paramPhiBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j)
    {
      for(k = 0; k < sizePhiBins; ++k)
      {
        fBinMap[paramEtaBins[j].GetDouble()].insert(paramPhiBins[k].GetDouble());
      }
    }
  }

  // for better performance we transform map of sets to parallel vectors:
  // vector< double > and vector< vector< double >* >
  for(itEtaBin = fBinMap.begin(); itEtaBin != fBinMap.end(); ++itEtaBin)
  {
    fEtaBins.push_back(itEtaBin->first);
    phiBins = new vector< double >(itEtaBin->second.size());
    fPhiBins.push_back(phiBins);
    phiBins->clear();
    for(itPhiBin = itEtaBin->second.begin(); itPhiBin != itEtaBin->second.end(); ++itPhiBin)
    {
      phiBins->push_back(*itPhiBin);
    }
  }

  // read energy fractions for different particles
  param = GetParam("EnergyFraction");
  size = param.GetSize();

  // set default energy fractions values
  fFractionMap.clear();
  fFractionMap[0] = make_pair(0.0, 1.0);

  for(i = 0; i < size/2; ++i)
  {
    paramFractions = param[i*2 + 1];

    ecalFraction = paramFractions[0].GetDouble();
    hcalFraction = paramFractions[1].GetDouble();

    fFractionMap[param[i*2].GetInt()] = make_pair(ecalFraction, hcalFraction);
  }

  // read min E value for timing measurement in ECAL
  fTimingEnergyMin = GetDouble("TimingEnergyMin",4.);
  // For timing
  // So far this flag needs to be false
  // Curved extrapolation not supported
  fElectronsFromTrack = false;

  // should be optimised depending calorimeter resolution
  fECalMinSignificance = GetDouble("ECalMinSignificance", 0.0);
  fHCalMinSignificance = GetDouble("HCalMinSignificance", 0.0);

  // switch on or off the dithering of the center of DualReadoutCalorimeter towers
  fSmearTowerCenter = GetBool("SmearTowerCenter", true);

  // switch on or off the log normal smearing (gaussian if false)
  fSmearLogNormal = GetBool("SmearLogNormal", true);

  // read resolution formulas
  fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output arrays
  fTowerOutputArray = ExportArray(GetString("TowerOutputArray", "towers"));
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "photons"));

  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "eflowTracks"));
  fEFlowPhotonOutputArray = ExportArray(GetString("EFlowPhotonOutputArray", "eflowPhotons"));
  fEFlowNeutralHadronOutputArray = ExportArray(GetString("EFlowNeutralHadronOutputArray", "eflowNeutralHadrons"));
}

//------------------------------------------------------------------------------

void DualReadoutCalorimeter::Finish()
{
  vector< vector< Double_t >* >::iterator itPhiBin;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  for(itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin)
  {
    delete *itPhiBin;
  }
}

//------------------------------------------------------------------------------

void DualReadoutCalorimeter::Process()
{
  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t ecalFraction, hcalFraction;
  Double_t ecalEnergy, hcalEnergy;
  Double_t sigma;
  Double_t energyGuess, energy;
  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector< Double_t >::iterator itEtaBin;
  vector< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  vector< Long64_t >::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fECalTowerFractions.clear();
  fHCalTowerFractions.clear();
  fECalTrackFractions.clear();
  fHCalTrackFractions.clear();

  // loop over all particles
  fItParticleInputArray->Reset();
  number = -1;
  fTowerRmax=0.;

  //cout<<"--------- new event ---------- "<<endl;

  while((particle = static_cast<Candidate*>(fItParticleInputArray->Next())))
  {
    const TLorentzVector &particlePosition = particle->Position;
    ++number;

    // compute maximum radius (needed in FinalizeTower to assess whether barrel or endcap tower)
    if (particlePosition.Perp() > fTowerRmax)
      fTowerRmax=particlePosition.Perp();

    pdgCode = TMath::Abs(particle->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
    {
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fECalTowerFractions.push_back(ecalFraction);
    fHCalTowerFractions.push_back(hcalFraction);

    if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), particlePosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // loop over all tracks
  fItTrackInputArray->Reset();
  number = -1;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trackPosition = track->Position;
    ++number;

    pdgCode = TMath::Abs(track->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
    {
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fECalTrackFractions.push_back(ecalFraction);
    fHCalTrackFractions.push_back(hcalFraction);

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), trackPosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  towerEtaPhi = 0;
  fTower = 0;
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits)
  {
    towerHit = (*itTowerHits);
    flags = (towerHit >> 24) & 0x00000000000000FFLL;
    number = (towerHit) & 0x0000000000FFFFFFLL;
    hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // finalize previous tower
      FinalizeTower();

      // create new tower
      fTower = factory->NewCandidate();

      phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

      // phi bins for given eta bin
      phiBins = fPhiBins[etaBin];

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5*(fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5*((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);

      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      fECalTowerEnergy = 0.0;
      fHCalTowerEnergy = 0.0;

      fECalTrackEnergy = 0.0;
      fHCalTrackEnergy = 0.0;
      fTrackEnergy = 0.0;

      fECalTrackSigma = 0.0;
      fHCalTrackSigma = 0.0;
      fTrackSigma = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fTowerTime = 0.0;
      fTowerTimeWeight = 0.0;

      fECalTowerTrackArray->Clear();
      fHCalTowerTrackArray->Clear();
      fTowerTrackArray->Clear();

    }

    // check for track hits
    if(flags & 1)
    {
      ++fTowerTrackHits;

      track = static_cast<Candidate*>(fTrackInputArray->At(number));
      momentum = track->Momentum;
      position = track->Position;

      ecalEnergy = momentum.E() * fECalTrackFractions[number];
      hcalEnergy = momentum.E() * fHCalTrackFractions[number];
      energy = ecalEnergy + hcalEnergy;

      if(ecalEnergy > fTimingEnergyMin && fTower)
      {
        if(fElectronsFromTrack)
        {
          fTower->ECalEnergyTimePairs.push_back(make_pair<Float_t, Float_t>(ecalEnergy, track->Position.T()));
        }
      }

      // in Dual Readout we do not care if tracks are ECAL of HCAL
      if(fECalTrackFractions[number] > 1.0E-9 || fHCalTrackFractions[number] > 1.0E-9)
      {
        fTrackEnergy += energy;
        // this sigma will be used to determine whether neutral excess is significant. We choose the resolution according to bthe higest deposited fraction (in practice had for charged hadrons and em for electrons)
        sigma = 0.0;
        if(fHCalTrackFractions[number] > 0)
          sigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());
        else
          sigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, momentum.E());

        if(sigma/momentum.E() < track->TrackResolution)
          energyGuess = ecalEnergy + hcalEnergy;
        else
          energyGuess = momentum.E();

        fTrackSigma += (track->TrackResolution)*energyGuess*(track->TrackResolution)*energyGuess;

        //add-----------------------------------------------------
        track->Eem = ecalEnergy;
        track->Ehad = hcalEnergy;
        //----------------------------------------------------------

        fTowerTrackArray->Add(track);
<<<<<<< HEAD


=======
>>>>>>> 7b7cc48d2c23d91a4377bbe83953fb58ae7aa317
      }
      else
      {
        fEFlowTrackOutputArray->Add(track);
      }

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;
    position = particle->Position;


    // fill current tower
    ecalEnergy = momentum.E() * fECalTowerFractions[number];
    hcalEnergy = momentum.E() * fHCalTowerFractions[number];

    fECalTowerEnergy += ecalEnergy;
    fHCalTowerEnergy += hcalEnergy;

    // assume combined timing measurements in ECAL/HCAL sections
    fTowerTime += (ecalEnergy + hcalEnergy) * position.T(); //sigma_t ~ 1/sqrt(E)
    fTowerTimeWeight += ecalEnergy + hcalEnergy;

    fTower->AddCandidate(particle);
    fTower->Position = position;
  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void DualReadoutCalorimeter::FinalizeTower()
{

  Candidate *track, *tower, *mother, *candidate;
  Double_t energy, pt, eta, phi, r, time;
  Double_t neutralEnergy;

  Double_t caloSigma, trackCaloSigma;
  Double_t neutralSignificance;
  Double_t neutralMinPFSignificance;

  Double_t weightTrack, weightCalo, bestEnergyEstimate, rescaleFactor;
  Bool_t isPureEM = false;

  TLorentzVector momentum;
  TFractionMap::iterator itFractionMap;

  Bool_t debug = false;
  if(!fTower) return;
  if (debug) cout<<"-----------------------------------------------------------------------"<<endl;
  if (debug) cout<<"New Tower: "<<fECalTowerEnergy<<","<<fHCalTowerEnergy<<","<<fHCalTowerEnergy<<","<<fTowerEta<<endl;

  if (debug) cout<<"   gen particles in tower :"<<fTower->GetCandidates()->GetEntriesFast()<<endl;
  TIter it1(fTower->GetCandidates());
  it1.Reset();
  while((candidate = static_cast<Candidate *>(it1.Next())))
  {
     //cout<<": " << <<endl;
     TLorentzVector mom = candidate->Momentum;
     if (debug) cout<<"      gen particle: " << candidate->PID<<"," << mom.E()<<"," << mom.Eta()<<"," <<mom.Phi()<<endl;
  }

  // if no hadronic energy, use ECAL resolution
  if (fHCalTowerEnergy <= 0)
  {
    energy = fECalTowerEnergy;
    caloSigma  = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    isPureEM = true;
    if (debug) cout<<"   using ECAL energy: "<<energy<<", "<<caloSigma<<endl;
  }

  // if hadronic fraction > 0, use HCAL resolution
  else
  {
    energy = fECalTowerEnergy + fHCalTowerEnergy;
    caloSigma  = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    if (debug) cout<<"   using HCAL energy: "<<energy<<", "<<caloSigma<<endl;
  }

  if (fSmearLogNormal)
    energy = LogNormal(energy, caloSigma);
  else
    //energy = TruncatedGaussian(energy, caloSigma);
    energy = gRandom->Gaus(energy, caloSigma);

  if (debug) cout<<"   smeared energy: "<<energy<<endl;

  if (energy < 0.) energy = 0.;
  // set tower energy to 0 when energy deposit is not significant

  if (isPureEM)
  {
    // estimate resolution from the measurement this time
    caloSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    energy = (energy > fECalMinSignificance * caloSigma) ? energy : 0.;
  }
  else
  {
    // estimate resolution from the measurement this time
    caloSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, energy);
    energy = (energy > fHCalMinSignificance * caloSigma) ? energy : 0.;
  }

  if (debug) cout<<"   smeared energy: "<<energy<<endl;

  // ---------------------------------------------------------------------------
  // compute calo tower properties
  // ---------------------------------------------------------------------------

  time = (fTowerTimeWeight < 1.0E-09) ? 0.0 : fTowerTime / fTowerTimeWeight;

  if(fSmearTowerCenter)
  {
    eta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
    phi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);
  }
  else
  {
    eta = fTowerEta;
    phi = fTowerPhi;
  }

  // check whether barrel or endcap tower

  // endcap
  if (TMath::Abs(fTower->Position.Pt() - fTowerRmax) > 1.e-06 && TMath::Abs(eta) > 0.){
    r = fTower->Position.Z()/TMath::SinH(eta);
  }
  // barrel
  else {
    r = fTower->Position.Pt();
  }

  time = (fTowerTimeWeight < 1.0E-09) ? 0.0 : fTowerTime / fTowerTimeWeight;


  fTower->Position.SetPtEtaPhiE(r, eta, phi, time);
  fTower->L = fTower->Position.Vect().Mag();
  fTower->Etrk = fTrackEnergy;

  // these are stored for debugging purposes, should not be used since they are
  // based on MC truth
  fTower->Eem = fECalTowerEnergy;
  fTower->Ehad = fHCalTowerEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  if (isPureEM)
  {
    // assume massless photon hypothesis
    fTower->PID = 22;
    pt = energy / TMath::CosH(eta);
    fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  }
  // if hadronic fraction > 0, use HCAL resolution
  else
  {

    // assume pion hypothesis for hadronic deposit. This can be corrected later by accessing particle energy in the output
    fTower->PID = 211;
    Double_t mass = 0.13957;
    Double_t p = (energy > mass) ? TMath::Sqrt(energy*energy - mass*mass) : 0.;
    pt = p / TMath::CosH(eta);
    fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  }

  if(energy > 0.0)
  {
    if(fTowerPhotonHits > 0 && fTowerTrackHits == 0)
    {
      fPhotonOutputArray->Add(fTower);
    }

    if (debug) cout<<"   creating tower with energy: "<<energy<<endl;
    if (debug) cout<<"   creating tower with PID: "<<fTower->PID<<endl;
    if (debug) cout<<"   creating tower with track energy: "<<fTower->Etrk<<endl;

    fTowerOutputArray->Add(fTower);
  }

  // ---------------------------------------------------------------------------
  // now do particle-flow
  // ---------------------------------------------------------------------------

  fTrackSigma = TMath::Sqrt(fTrackSigma);
  neutralEnergy = max( (energy - fTrackEnergy) , 0.0);

  if (isPureEM)
  {
    neutralMinPFSignificance = fECalMinSignificance;
  }
  else
  {
    neutralMinPFSignificance = fHCalMinSignificance;
  }

  // combined track calo resolution
  trackCaloSigma = TMath::Sqrt(fTrackSigma*fTrackSigma + caloSigma*caloSigma);
  neutralSignificance = (trackCaloSigma > 0) ? neutralEnergy / trackCaloSigma : 0.;

  if (debug) cout<<"Doing PF here: "<<endl;
  if (debug) cout<<"   track energy: "<<fTrackEnergy<<endl;
  if (debug) cout<<"   calo energy: "<<energy<<endl;
  if (debug) cout<<"   neutral energy: "<<neutralEnergy<<endl;
  if (debug) cout<<"   track sigma: "<<fTrackSigma<<endl;
  if (debug) cout<<"   calo sigma: "<<caloSigma<<endl;
  if (debug) cout<<"   track calo sigma: "<<trackCaloSigma<<endl;
  if (debug) cout<<"   neutral significance: "<<neutralSignificance<<endl;


  // now do case where at least one track points to tower and the nuetral excess is signficant
  // i.e pi+ and neutron or electron and photon hitting same tower
  if ( neutralSignificance > neutralMinPFSignificance )
  {

    tower = static_cast<Candidate*>(fTower->Clone());
    if (isPureEM)
    {
      tower->Eem = neutralEnergy;
      tower->Ehad = 0.0;
      tower->PID = 22;
      pt = neutralEnergy / TMath::CosH(eta);
      tower->Momentum.SetPtEtaPhiE(pt, eta, phi, neutralEnergy);
      fEFlowPhotonOutputArray->Add(tower);
    }
    else
    {
      tower->Eem = 0;
      tower->Ehad = neutralEnergy;
      tower->PID = 130;
      Double_t mass = 0.497611;
      Double_t p = (neutralEnergy > mass) ? TMath::Sqrt(neutralEnergy*neutralEnergy - mass*mass) : 0.;
      pt = p / TMath::CosH(eta);
      if (p > 0)
      {
        tower->Momentum.SetPtEtaPhiE(pt, eta, phi, neutralEnergy);
        fEFlowNeutralHadronOutputArray->Add(tower);
      }
    }

    if (debug) cout<<"       creating neutral excess with energy, eta, phi: "<<neutralEnergy<<","<<eta<<","<<phi<<endl;
    if (debug) cout<<"       creating neutral excess with PID: "<<tower->PID<<endl;
    if (debug) cout<<"       creating neutral excess with track energy: "<<tower->Etrk<<endl;
    if (debug) cout<<"       "<<endl;

    // now clone tracks
    fItTowerTrackArray->Reset();
    while((track = static_cast<Candidate*>(fItTowerTrackArray->Next())))
    {
      mother = track;
      track = static_cast<Candidate*>(track->Clone());
      track->AddCandidate(mother);
      
      //add-------------------------------------------------------------------
      Double_t Eema = track->Eem;
      Double_t Ehada = track->Ehad;
      Double_t EcalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, Eema);
      Double_t HcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, Ehada);
      track->Eem = LogNormal(Eema, EcalSigma);
      track->Ehad = LogNormal(Ehada, HcalSigma);
      //-------------------------------------------------------------------

      fEFlowTrackOutputArray->Add(track);
    }

  }

  // now do case where at track points to tower and the neutral excess is NOT signficant
  // if neutral excess is not significant, rescale eflow tracks, such that the total
  // charged equals the best measurement given by the DualReadoutCalorimeter and tracking
  else if (fTrackEnergy > 0)
  {
    if (debug) cout<<"   no significant neutral excess found:"<<endl;
    if (debug) cout<<"   neutral energy: "<<energy<<", "<<fTrackEnergy<<", "<<neutralEnergy<<endl;
    if (debug) cout<<"       "<<endl;

    weightTrack = (fTrackSigma > 0.0) ? 1 / (fTrackSigma*fTrackSigma) : 0.0;
    weightCalo  = (caloSigma > 0.0) ? 1 / (caloSigma*caloSigma) : 0.0;

    bestEnergyEstimate = (weightTrack*fTrackEnergy + weightCalo*energy) / (weightTrack + weightCalo);
    rescaleFactor = bestEnergyEstimate/fTrackEnergy;

    //rescale tracks
    fItTowerTrackArray->Reset();
    while((track = static_cast<Candidate*>(fItTowerTrackArray->Next())))
    {
      mother = track;
      track = static_cast<Candidate *>(track->Clone());
      track->AddCandidate(mother);

      //add-------------------------------------------------
      Double_t ecalf = track->Eem/track->Momentum.E();
      Double_t hcalf = track->Ehad/track->Momentum.E();
      //----------------------------------------------------------

      track->Momentum.SetPtEtaPhiM(track->Momentum.Pt()*rescaleFactor, track->Momentum.Eta(), track->Momentum.Phi(), track->Momentum.M());
<<<<<<< HEAD


      //add-------------------------------------------------------------
      Double_t Eema = ecalf*track->Momentum.E();
      Double_t Ehada = hcalf*track->Momentum.E();
      Double_t EcalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, Eema);
      Double_t HcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, Ehada);
      track->Eem = LogNormal(Eema, EcalSigma);
      track->Ehad = LogNormal(Ehada, HcalSigma);
      //---------------------------------------------------------------
=======
      if (debug) cout<<"  track Momentum: "<< track->PID<<", "<<track->Momentum.Pt()<<", "<<track->Momentum.Eta()<<", "<<track->Momentum.M()<<endl;
>>>>>>> 7b7cc48d2c23d91a4377bbe83953fb58ae7aa317
      fEFlowTrackOutputArray->Add(track);
    }
  }

}


//------------------------------------------------------------------------------

Double_t DualReadoutCalorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
    a = TMath::Log(mean) - 0.5*b*b;

    return TMath::Exp(a + b*gRandom->Gaus(0.0, 1.0));
  }
  else
  {
    return 0.0;
  }
}

//------------------------------------------------------------------------------

Double_t DualReadoutCalorimeter::TruncatedGaussian(Double_t mean, Double_t sigma)
{
  Double_t result = -1;
  if(mean > 0.0)
  {
    while (result < 0.0)
    {
      result = gRandom->Gaus(mean, sigma);
    }
    return result;
  }
  else
  {
    return 0.0;
  }
}
