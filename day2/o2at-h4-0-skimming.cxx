// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task contains the individual steps that are to be taken
///        in the first part of the tutorial. These are 6 steps, and at the end,
///        the participant is expected to have a pT resolution histogram
///        as a function of transverse momentum.
/// \author
//	Johanna
//	Additional tasks required:
//	– o2-analysis-timestamp
//	– o2-analysis-event-selection
//	– o2-analysis-multiplicity-table
//	– o2-analysis-track-propagation
//	– o2-analysis-trackselection
//	– o2-analysis-hf-track-index-skims-creator
//	– o2-analysis-hf-candidate-creator-2prong
//
//to run:
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 1
//<- starting point, define the derived table to be stored
// this can be done in a separated header file, but for semplicity we do it in
// the same file here
// with the only table definition, there is nothing to be executed

namespace o2::aod
{
namespace mytable
{
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);           //!
DECLARE_SOA_COLUMN(InvMassD0bar, invMassD0bar, float);     //!
DECLARE_SOA_COLUMN(Pt, pt, float);                         //!
DECLARE_SOA_COLUMN(CosinePointing, cosinePointing, float); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);            //!
} // namespace mytable

DECLARE_SOA_TABLE(MyTable, "AOD", "MYTABLE", //!
                  mytable::InvMassD0,
                  mytable::InvMassD0bar,
                  mytable::Pt,
                  mytable::CosinePointing,
                  mytable::CollisionId)

} // namespace o2::aod
 
struct ReadHFCandidates{// simple workflow that loops over HF 2-prong candidates 
  
  void process(aod::HfCandProng2 const& cand2Prongs){
    //loop over HF 
    for(auto& cand : cand2Prongs){
      //check first of HF 2-prong candidate is tagged as D0 -> this will decay into pion and kaon ? or pion or kaon ?
      bool isD0Sel = TESTBIT(cand.hfflag(), aod::hf_cand_prong2::DecayType::D0ToPiK);

      if(!isD0Sel){
        continue;
      }

      auto invMassD0 = InvMassD0(cand);
      auto invMassD0bar = InvMassD0bar(cand);

      LOG(debug) << "Candidate with mass(D0) = " << invMassD0
	         << ", mass(D0bar) = " << invMassD0bar
		 << ", pt = " << cand.pt()
		 << ", cos(theta_P) = " << cand.cpa();
    }

  }

};

//Step 2 - Read in HF candidates-nothing will happen to the produced table -> create a task in next step
struct ProduceDerivedTable{// workflow that loops over HF 2-prong candidates and fills a derived table

  Produces<aod::MyTable> tableWithDzeroCandidates;//above is the declaration of the table

  void process(aod::HfCandProng2 const& cand2Prongs, aod::Tracks const&){
  
    //loop over HF 
    for(auto& cand : cand2Prongs){
      //check first of HF 2-prong candidate is tagged as D0 -> this will decay into pion and kaon ? or pion or kaon ?
      bool isD0Sel = TESTBIT(cand.hfflag(), aod::hf_cand_prong2::DecayType::D0ToPiK);

      //select only D0 candidates with pT > 4 GeV/c
      if(!isD0Sel || cand.pt() < 4.){
        continue;
      }

      auto invMassD0 = InvMassD0(cand);
      auto invMassD0bar = InvMassD0bar(cand);

      LOG(debug) << "Candidate with mass(D0) = " << invMassD0
                 << ", mass(D0bar) = " << invMassD0bar
                 << ", pt = " << cand.pt()
                 << ", cos(theta_P) = " << cand.cpa();
      //retrieve also the event index from one of the daughters
      auto dauTrack = cand.index0_as<aod::Tracks>(); //positive daughter
	
      tableWithDzeroCandidates(invMassD0, invMassD0bar, cand.pt(), cand.cpa(), dauTrack.collisionId());//is produced following the structure in the first lines.
    }

  }

};

//Step 3
struct ReadDerivedTable{//workflow that reads derived tabke and fills histograms

  HistogramRegistry registry{
    "registry",
     {
     	{"hMassD0",";#it{M}(K#pi) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{300, 1.75, 2.05}}}},
	{"hMassD0bar", ";#it{M}(#piK) (GeV/#it{c}^{2};counts)", {HistType::kTH1F, {{300, 1.75, 2.05}}}},
	{"hPt", ";#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{50,0,50}}}},
	{"hCosp", ";cos(#vartheta_{p}) ;counts", {HistType::kTH1F, {{100,0.8,1.}}}}
     }
  };

  void process(aod::MyTable const& cand2Prongs){
     //loop over 2-prong candidates
    for(auto& cand : cand2Prongs){
      registry.fill(HIST("hMassD0"), cand.invMassD0());
      registry.fill(HIST("hMassD0bar"), cand.invMassD0bar());
      registry.fill(HIST("hPt"), cand.pt());
      registry.fill(HIST("hCosp"), cand.cosinePointing());
    }
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
     // adaptAnalysisTask<ProduceDerivedTable>(cfgc), -> I produced that table with the OutputDirector.json that was in the initial run_*.sh
      adaptAnalysisTask<ReadDerivedTable>(cfgc)
  };
}

