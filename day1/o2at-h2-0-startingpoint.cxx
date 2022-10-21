// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task contains the individual steps that are to be taken
///        in the second part of the tutorial. These are 5 steps, and at the end,
///        the participant is expected to have a two-particle correlation spectrum.
/// \author
//	Johanna
//
//	to run: o2-analysistutorial-o2at-h3-0-startingpoint --aod-file AO2D-LHC21k6-2.root | o2-analysis-track-propagation | o2-analysis-timestamp -b
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//This is an example of a conveient declaration of "using"
//WARNING: THIS IS NEW
using MyCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

//Step 1-2
struct filterexample {
  Filter etaFilter = nabs(aod::track::eta) < 0.5f;
  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogram", "etaHistogram", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogram", "ptHistogram", {HistType::kTH1F, {{nBins, 0., 10.0}}}}
    }
  };

  void process(aod::Collision const& collision, soa::Filtered<MyCompleteTracks> const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    //This will take place once per event!
    for (auto& track : tracks) {
      //"Good old" imperative cuts
      if( fabs(track.eta()) > 0.5 ) continue;
      if( track.tpcNClsCrossedRows() < 70 ) continue;//can't filter on dynamic
      if( fabs(track.dcaXY()) > .2 ) continue;
      registry.get<TH1>(HIST("etaHistogram"))->Fill(track.eta());
      registry.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
    }
  }
};

//Step 3
struct partitionexample {
  Partition<MyCompleteTracks> leftTracks = aod::track::eta < 0.0f;
  Partition<MyCompleteTracks> rightTracks = aod::track::eta >= 0.0f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxZ", "hVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHistogramright", "etaHistogramright", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"etaHistogramleft", "etaHistogramleft", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"ptHistogramright", "ptHistogramright", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"ptHistogramleft", "ptHistogramleft", {HistType::kTH1F, {{nBins, 0., 10.0}}}},

    }
  };

  void process(aod::Collision const& collision,/* soa::Filtered<*/MyCompleteTracks/*>*/ const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());

    //partitions are not grouped by default
    auto leftTracksGrouped = leftTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto rightTracksGrouped = rightTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    for(auto& track : leftTracksGrouped){// only for this subset
        registry.get<TH1>(HIST("etaHistogramleft"))->Fill(track.eta());
	registry.get<TH1>(HIST("ptHistogramleft"))->Fill(track.pt());
    }
    for(auto& track : rightTracksGrouped){// only for this subset
        registry.get<TH1>(HIST("etaHistogramright"))->Fill(track.eta());
        registry.get<TH1>(HIST("ptHistogramright"))->Fill(track.pt());
    }
  }
};

//Step 4
struct partandfilterexample {
  //all defined filters are applied 
  Filter etaFilter = nabs(aod::track::eta) < 0.5f;
  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;
  using MyFilteredTracks = soa::Filtered<MyCompleteTracks>;

  Partition<MyFilteredTracks> leftTracks = aod::track::eta < 0.0f;//adjust according to filter
  Partition<MyFilteredTracks> rightTracks = aod::track::eta >= 0.0f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"HVtxZ", "HVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"etaHright", "etaHright", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"etaHleft", "etaHleft", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"ptHright", "ptHright", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"ptHleft", "ptHleft", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
    }
  };

  void process(aod::Collision const& collision, MyFilteredTracks const& tracks) //<- this is the main change
  {
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("HVtxZ"))->Fill(collision.posZ());

    //partitions are not grouped by default
    auto leftTracksGrouped = leftTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto rightTracksGrouped = rightTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    for(auto& track : leftTracksGrouped){// only for this subset
	if(track.tpcNClsCrossedRows() <70 ) continue;// can't filter on dynamic
        registry.get<TH1>(HIST("etaHleft"))->Fill(track.eta());
        registry.get<TH1>(HIST("ptHleft"))->Fill(track.pt());
    }
    for(auto& track : rightTracksGrouped){// only for this subset
        if(track.tpcNClsCrossedRows() <70 ) continue;//can't filter on dynamic
        registry.get<TH1>(HIST("etaHright"))->Fill(track.eta());
        registry.get<TH1>(HIST("ptHright"))->Fill(track.pt());
    }
  }
};

//Step 4 & 5 - ofc I could include QA and the secret step...but that is for a rainy day...
struct twoparticlecorrelation {
  //all defined filters are applied
  Filter trackFilter = nabs(aod::track::eta) < 0.8f && aod::track::pt > 2.0f;
  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;
  using MyFilteredTracks = soa::Filtered<MyCompleteTracks>;

  Partition<MyFilteredTracks> triggerTrack = aod::track::pt > 4.0f && aod::track::pt < 10.0f;//adjust according to filter
  Partition<MyFilteredTracks> assocTrack = aod::track::pt < 4.0f;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxZ", "hVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"correlationFunction", "correlationFunction", {HistType::kTH2F, {{nBins, -10., 10.}, {nBins, -10,10}}}},
      {"correlationFunctionComb", "correlationFunctionComb", {HistType::kTH2F, {{nBins, -10., 10.}, {nBins,-10,10}}}},
      {"correlationPhi", "correlationPhi", {HistType::kTH1F, {{nBins, -10., 10.}}}},
      {"correlationEta", "correlationEta", {HistType::kTH1F, {{nBins, -10., 10.}}}},
      {"etaHtrigger", "etaHtrigger", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"etaHassoc", "etaHassoc", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"phiHtrigger", "phiHtrigger", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"phiHassoc", "phiHassoc", {HistType::kTH1F, {{nBins, -1, 1}}}},
      {"ptHtrigger", "ptHtrigger", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"ptHassoc", "ptHassoc", {HistType::kTH1F, {{nBins, 0., 10.0}}}}
     
     }
  };

  Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) {
    //To be completely sure, use inner products to get delta phi = phi1-phi2
    Double_t x1, y1, x2, y2;
    x1 = TMath::Cos( phi1 );
    y1 = TMath::Sin( phi1 );
    x2 = TMath::Cos( phi2 );
    y2 = TMath::Sin( phi2 );
    Double_t lInnerProd = x1*x2 + y1*y2;
    Double_t lVectorProd = x1*y2 - x2*y1;
    Double_t lReturnVal = 0;
    if( lVectorProd > 1e-8 ) {
      lReturnVal = TMath::ACos(lInnerProd);
    }
    if( lVectorProd < -1e-8 ) {
      lReturnVal = -TMath::ACos(lInnerProd);
    }
    if( lReturnVal < -TMath::Pi()/2. ) {
      lReturnVal += 2.*TMath::Pi();
    }
    return lReturnVal;
  } 

  void process(aod::Collision const& collision, MyFilteredTracks const& tracks) //<- this is the main change and maybe my problem .. or the internet here...
  {
    //Fill the event counter
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
    auto triggerTracks = triggerTrack->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto assocTracks = assocTrack->sliceByCached(aod::track::collisionId, collision.globalIndex());
    //two-particle correlation - still manually
    for(auto& trackTrigger : triggerTracks){// only for this subset
      if(trackTrigger.tpcNClsCrossedRows() <70 ) continue;// can't filter on dynamic
      for(auto& trackAssoc : assocTracks){// only for this subset
        if(trackAssoc.tpcNClsCrossedRows() <70 ) continue;//can't filter on dynamic
          registry.get<TH2>(HIST("correlationFunction"))->Fill(ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ), trackTrigger.eta() - trackAssoc.eta() );
        }
    }
    
    //two-particle correlation - via combinations - fancy - but i cannot see the 
    for(auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(triggerTracks, assocTracks)) ){
      if(trackTrigger.tpcNClsCrossedRows() <70 ) continue;// can't filter on dynamic
      if(trackAssoc.tpcNClsCrossedRows() <70 ) continue;//can't filter on dynamic
      registry.get<TH2>(HIST("correlationFunctionComb"))->Fill(ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ), trackTrigger.eta() - trackAssoc.eta() );
      registry.get<TH1>(HIST("correlationPhi"))->Fill(ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() )); 
      registry.get<TH1>(HIST("correlationEta"))->Fill(trackTrigger.eta() - trackAssoc.eta() );
      //control hists
      registry.get<TH1>(HIST("etaHtrigger"))->Fill(trackTrigger.eta() );
      registry.get<TH1>(HIST("etaHassoc"))->Fill(trackAssoc.eta() );
      registry.get<TH1>(HIST("phiHtrigger"))->Fill(trackTrigger.phi() );
      registry.get<TH1>(HIST("phiHassoc"))->Fill(trackAssoc.phi() );
      registry.get<TH1>(HIST("ptHtrigger"))->Fill(trackTrigger.pt() );
      registry.get<TH1>(HIST("ptHassoc"))->Fill(trackAssoc.pt() );

    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<filterexample>(cfgc),
    adaptAnalysisTask<partitionexample>(cfgc),
    adaptAnalysisTask<partandfilterexample>(cfgc),
    adaptAnalysisTask<twoparticlecorrelation>(cfgc)
  };
}
