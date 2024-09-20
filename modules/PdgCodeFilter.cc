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

/** \class PdgCodeFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/PdgCodeFilter.h"

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

PdgCodeFilter::PdgCodeFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PdgCodeFilter::~PdgCodeFilter()
{
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Init()
{

  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  fEnMin = GetDouble("EnMin", 0.0);

  fMassRes = GetDouble("MassRes", -1.0);

  fNP      = GetDouble("NP", -1);

  fInvert = GetBool("Invert", false);

  fRequireStatus = GetBool("RequireStatus", false);
  fStatus = GetInt("Status", 1);

  fRequireCharge = GetBool("RequireCharge", false);
  fCharge = GetInt("Charge", 1);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  param = GetParam("PdgCode");
  size = param.GetSize();

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  fOutputArray1 = ExportArray(GetString("OutputArray1", "WoMuonPair"));
  fOutputArray2 = ExportArray(GetString("OutputArray2", "MuonPair"  ));
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Process()
{
  Candidate *candidate;
  Candidate *candidate1, *can1=NULL;
  Candidate *candidate2, *can2=NULL;
  Bool_t pass;
  Int_t pdgCode;
  Double_t pt, en;
  Int_t pdgCode1, pdgCode2;
  Double_t pt1, en1, pt2, en2;
  fItInputArray->Reset();
  if ( fNP==2 && fMassRes>0) {
	  //cout<< "I am here "<<fMassRes<<" "<<fNP <<" "<<fEnMin<<endl;
	  Double_t dm = 999999;
	  while((candidate1 = static_cast<Candidate*>(fItInputArray->Next()))){
		  pdgCode1 = candidate1->PID;
		  const TLorentzVector &candidateMomentum1 = candidate1->Momentum;
		  en1 = candidateMomentum1.E();

		  if(fRequireStatus && (candidate1->Status != fStatus)) continue;
		  if(fRequireCharge && (candidate1->Charge != fCharge)) continue;
		  if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode1) == fPdgCodes.end() || en1<fEnMin) continue;

		  while((candidate2 = static_cast<Candidate*>(fItInputArray->Next()))){
			  if( candidate2 == candidate1) continue;
			  pdgCode2 = candidate2->PID;
			  const TLorentzVector &candidateMomentum2 = candidate2->Momentum;
			  en2 = candidateMomentum2.E();

			  if(fRequireStatus && (candidate2->Status != fStatus)) continue;
			  if(fRequireCharge && (candidate2->Charge != fCharge)) continue;

			  if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode2) == fPdgCodes.end() || en2<fEnMin) continue;

			  double mass  = fabs((candidateMomentum1 + candidateMomentum2).M()-fMassRes);
			  if( mass<dm  ){
				  can1 = candidate1;
				  can2 = candidate2;
				  dm   = mass;
			  }
		  }
	  }

	  if(can1 && can2) {
		  fOutputArray2->Add(can1);
		  fOutputArray2->Add(can2);
	  }

	  fItInputArray->Reset();
	  TIterator *fItOutputArray = fOutputArray2->MakeIterator();
	  Int_t number=0;
	  while((candidate1 = static_cast<Candidate*>(fItInputArray->Next()))){
		  bool saved = false;
		  fItOutputArray->Reset();
		  while((candidate2 = static_cast<Candidate*>(fItOutputArray->Next())) ){
			  if ( candidate2 == candidate1 ) saved = true;
		  }
		  if ( ! saved ) { 
			  fOutputArray1->Add(candidate1); 
			  //cout<<"saved "<<number++<<endl;
		  }
	  }
  }else{
	  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
	  {
		  pdgCode = candidate->PID;
		  const TLorentzVector &candidateMomentum = candidate->Momentum;
		  pt = candidateMomentum.Pt();
		  en = candidateMomentum.E();

		  if(pt < fPTMin) continue;
		  if(fRequireStatus && (candidate->Status != fStatus)) continue;
		  if(fRequireCharge && (candidate->Charge != fCharge)) continue;

		  pass = kTRUE;
		  if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) != fPdgCodes.end() ) pass = kFALSE;
		  if(fInvert) pass = !pass;
		  if(pass) fOutputArray1->Add(candidate);
		  else     fOutputArray2->Add(candidate);
	  }
  }
}

