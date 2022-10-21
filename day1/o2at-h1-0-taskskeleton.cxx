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
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author
///	Johanna
//	
//	to run: o2-analysistutorial-o2at-h1-0-taskskeleton --aod-file AO2D-LHC21k6-2.root | o2-analysis-track-propagation | o2-analysis-timestamp -b
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;

//SECRET STEP: missing but planned as gap activity for the journal club prep.
//This is an empty analysis skeleton: was - now we have 1-7 from pT-resolution tutorial completed! 
struct taskskeleton {
  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -1., +1)};
  OutputObj<TH1F> ptHistogram{TH1F("ptHistogram","ptHistogram", 100, 0., 10.0)};

  //Configurable number of bins
  Configurable<int> nBinsEta{"nBinsEta",100, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt",100, "N bins in pT histo"};
  Configurable<int> nBinsZ{"nBinsZ",100,"N bins in z vertex histo"};

  //incomplete histogram definitions
  OutputObj<TH1F> etaHist{"etaHist"};
  OutputObj<TH1F> ptHist{"ptHist"};

  void init(InitContext const&)//example for different histogram initialisation 
  {//to complete definitioan of OutputObj
   etaHist.setObject(new TH1F("etaHist","etaHist", nBinsEta, -1,1));
   ptHist.setObject(new TH1F("ptHist","ptHist",nBinsPt,0.,10));
  }			      

  HistogramRegistry registry{//example of histogram registry -> just to understand what o2 can do for u
     "registry",
      {
	 {"hVertexZ", "hVertexZ", {HistType::kTH1F,{{nBinsZ,-10,10}}}},
	 {"hResolution", "hResolution", {HistType::kTH2F, {{nBinsPt,0,10},{nBinsZ,-5,5}}}},
	 {"etaH_TPC", "etaH_TPC", {HistType::kTH1F,{{nBinsEta,-1,1}}}},//adding _TPC to indicate TPC specific selection in process function
	 {"ptH_TPC", "ptH_TPC", {HistType::kTH1F, {{nBinsPt,0,10}}}},
         {"etaH_PV", "etaH_PV", {HistType::kTH1F,{{nBinsEta,-1,1}}}},//adding _PV to indicate additional Primary Vertex specific selection in process function
         {"ptH_PV", "ptH_PV", {HistType::kTH1F, {{nBinsPt,0,10}}}}
      }
  };
	
  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::McParticles const&)
  {//https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html to find what the tables can do for u
     //Fill event counter
     //This happens once per event !
     registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
     for (auto& track : tracks) {
	etaHistogram->Fill(track.eta());
	ptHistogram->Fill(track.pt());
	etaHist->Fill(track.eta());
	ptHist->Fill(track.pt());
	
	if( track.tpcNClsCrossedRows() < 70) continue;// TPC selection = skip messy TPC tracks
	registry.get<TH1>(HIST("etaH_TPC"))->Fill(track.eta());
	registry.get<TH1>(HIST("ptH_TPC"))->Fill(track.pt());
	if( fabs(track.dcaXY()) > 0.2) continue;// PV selection = kepp only tracks that point to the primary vertex
	registry.get<TH1>(HIST("etaH_PV"))->Fill(track.eta());
        registry.get<TH1>(HIST("ptH_PV"))->Fill(track.pt()); 
	//Resolve MC track 
	auto mcParticle = track.mcParticle_as<aod::McParticles>();
	float delta = track.pt() - mcParticle.pt();
	registry.get<TH2>(HIST("hResolution"))->Fill(track.pt(),delta);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskskeleton>(cfgc)
  };
}
