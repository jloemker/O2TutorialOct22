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
/// \brief this is a starting point for the third session of the tutorial
/// \author
//	Johanna
//
//	to run (at least Step 1,2):o2-analysistutorial-o2at-h3-0-startingpoint --configuration json://config.json | o2-analysis-track-propagation --configuration json://config.json | o2-analysis-timestamp --configuration json://config.json | o2-analysis-lf-lambdakzerobuilder --configuration json://config.json | o2-analysis-event-selection --configuration json://config.json -b
//
//	to run (step 3): o2-analysis-pid-tpc | o2-analysis-multiplicity-table | o2-analysistutorial-o2at-h3-0-startingpoint --configuration json://config.json | o2-analysis-track-propagation --configuration json://config.json | o2-analysis-timestamp --configuration json://config.json | o2-analysis-lf-lambdakzerobuilder --configuration json://config.json | o2-analysis-event-selection --configuration json://config.json -b
//	to run (step 5):  o2-analysis-pid-tpc --configuration json://config.json | o2-analysis-multiplicity-table --configuration json://config.json | o2-analysistutorial-o2at-h3-0-startingpoint --configuration json://config.json | o2-analysis-track-propagation --configuration json://config.json | o2-analysis-timestamp --configuration json://config.json | o2-analysis-lf-lambdakzerobuilder --configuration json://config.json | o2-analysis-event-selection --configuration json://config.json -b
//
//\since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;
using MyTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 2 data stores Tracks in AO2D
using MyTracksRun3 = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 3 data stores TracksIU (inner most update) in AO2D
using LabeledV0s = soa::Join<aod::V0Datas, aod::McV0Labels>;

//STEP 1
struct vzeroexample {

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}}//why these exact numbers ?
    }
  };

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& V0s) //<- this is the main change
  {
    //Basic event selection (all helper tasks are now included!)
    //https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
    //o2::aod::evsel::Sel7 		sel7 	bool 	Event selection decision based on V0A & V0C
    //o2::aod::evsel::Sel8 		sel8 	bool 	Event selection decision based on TVX)
    // (...)
    if(!collision.sel8()){
	    return;
    }
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    for(auto& v0 : V0s){
      registry.fill(HIST("hMassK0Short"), v0.mK0Short());
    }
  }
};


//STEP 3
//Now adding particle identification by using the tpcNSigma track variable
struct vzerofilterexample {
  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  //Selection criteria: 5 basic V0 selection criteria!
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA" };// double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", -1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", -1, "DCA Pos To PV"};//why would they also put -1 here ? 
  Configurable<float> v0radius{"v0radius", 0.5, "v0radius"};

  //Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daus only
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&
	  	       nabs(aod::v0data::dcanegtopv) > dcanegtopv&&
		       aod::v0data::dcaV0daughters < dcav0dau;
  //Define tracks for two separate process functions - can also be put outside of this scope !
 // using MyTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 2 data stores Tracks in AO2D
 // using MyTracksRun3 = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 3 data stores TracksIU (inner most update) in AO2D
  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxZ", "hVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hMK0Short", "hMK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},//why these exact numbers - around mass peak i assume :P
      {"hMLambda", "hMLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMAntiLambda", "hMAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}}
    }
  };

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& V0s, MyTracks const& tracks)
  {
    //Basic event selection
    // (...)
    // Now we need additionally o2-analysis-pid-tpc | o2-analysis-multiplicity-table
    if(!collision.sel8()){
            return;
    }
    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
    for(auto& v0 : V0s){
      float nsigma_pos_proton = TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPr());
      float nsigma_neg_proton = TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPr());
      float nsigma_pos_pion = TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      float nsigma_neg_pion = TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi());

      if(v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa ){
        if( nsigma_pos_pion < 4 && nsigma_neg_pion < 4 ){
	  registry.fill(HIST("hMK0Short"), v0.mK0Short());
	}
	if( nsigma_pos_proton < 4 && nsigma_neg_proton < 4){
          registry.fill(HIST("hMLambda"), v0.mLambda());
	}
	if( nsigma_pos_pion < 4 && nsigma_neg_proton < 4){
          registry.fill(HIST("hMAntiLambda"), v0.mAntiLambda());	
	} 
      }
    }
  }
};

struct vzerotemplateexample {
  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  //Selection criteria: 5 basic V0 selection criteria!
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA" };// double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", -1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", -1, "DCA Pos To PV"};//why would they also put -1 here ? 
  Configurable<float> v0radius{"v0radius", 0.5, "v0radius"};

  //Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daus only
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&
                       nabs(aod::v0data::dcanegtopv) > dcanegtopv&&
                       aod::v0data::dcaV0daughters < dcav0dau;
  //Define tracks for two separate process functions - can ofcourse happen outside the struct...
 /* using MyTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 2 data stores Tracks in AO2D
  using MyTracksRun3 = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 3 data stores TracksIU (inner most update) in AO2D
  using LabeledV0s = soa::Join<aod::V0Datas, aod::McV0Labels>;
 */ // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxZ", "hVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hMK0Short", "hMK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},//why these exact numbers - around mass peak i assume :P
      {"hMLambda", "hMLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMAntiLambda", "hMAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMassTrueK0Short", "hMassTrueK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},
      {"hMassTrueLambda", "hMassTrueLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMassTrueAntiLambda", "hMassTrueAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}}
    }
  };

  //to define the process module in a generic way-but this spoils all other implementations :(
  template <class TMyTracks, typename TV0>
  void processV0Candidate(TV0 const& v0, float const& pvx, float const& pvy, float const& pvz)
  {//function to process a vzero candidate freely, actually with right track type !
   auto posTrackCast = v0.template posTrack_as<TMyTracks>();
   auto negTrackCast = v0.template negTrack_as<TMyTracks>();

   float nsigma_pos_proton = TMath::Abs(posTrackCast.tpcNSigmaPr());
   float nsigma_neg_proton = TMath::Abs(posTrackCast.tpcNSigmaPr());//shouldn't it be negTrackCast ??
   float nsigma_pos_pion = TMath::Abs(negTrackCast.tpcNSigmaPi());//shouldn't it be posTrackCast ??
   float nsigma_neg_pion = TMath::Abs(negTrackCast.tpcNSigmaPi());
   
   if( v0.v0radius() > v0radius && v0.v0cosPA(pvx, pvy, pvz) > v0cospa ){
     if( nsigma_pos_pion < 4 && nsigma_neg_pion < 4){
       registry.fill(HIST("hMK0Short"), v0.mK0Short());
     }
     if( nsigma_pos_proton < 4 && nsigma_neg_proton < 4){
       registry.fill(HIST("hMLambda"), v0.mLambda());
     }
     if( nsigma_pos_pion < 4 && nsigma_neg_proton < 4){
       registry.fill(HIST("hMAntiLambda"), v0.mAntiLambda());
     }   
     // for this i need an additional .h file and another subscription in my pipeline !
     //An adding the final STEP *5* - lambdakzerobuilder is able to do MC association for me #wow !
     if( v0.has_mcParticle()){//soem association was made !
       auto v0mcparticle = v0.mcParticle();
       //Check particle PDG code to see if this is the one you want
       if( v0mcparticle.pdgCode() == 310 ) registry.fill(HIST("hMassTrueK0Short"), v0.mK0Short());
       if( v0mcparticle.pdgCode() == 3122 ) registry.fill(HIST("hMassTrueLambda"), v0.mLambda());// wtf why the same odg code for lamda and anti lambda ??
       if( v0mcparticle.pdgCode() == 3122 ) registry.fill(HIST("hMassTrueAntiLambda"), v0.mAntiLambda());
     }
   }
  }

  //Process function for old data (run2)
  void processRun2(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<LabeledV0s> const& V0s, MyTracksRun2 const& tracks, aod::McParticles const&)
  {
    if( !collision.sel7()){//sel7=bool -> Event selection decision based on V0A & V0C + all helper tasks required for this !
      return;
    }
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
    for( auto& v0 : V0s){
      processV0Candidate<MyTracksRun2>(v0, collision.posX(), collision.posY(), collision.posZ() );
    }
  }
  PROCESS_SWITCH(vzerotemplateexample, processRun2, "Process Run 2 data", false);
  //Define process function for run3 data
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<LabeledV0s> const& V0s, MyTracksRun3 const& tracks, aod::McParticles const&)
  {
    if( !collision.sel8()){//sel8=bool -> Event selection decision based on TVX) ! Event selection is based on a different detector part for run2 / run3 !
      return;
    }
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
    for( auto& v0 : V0s){
      processV0Candidate<MyTracksRun3>(v0, collision.posX(), collision.posY(), collision.posZ());
    }
  }
  PROCESS_SWITCH(vzerotemplateexample, processRun3, "Process Run 3 data", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
   // adaptAnalysisTask<vzeroexample>(cfgc)
   // adaptAnalysisTask<vzerofilterexample>(cfgc)
    adaptAnalysisTask<vzerotemplateexample>(cfgc)    
  };
}
