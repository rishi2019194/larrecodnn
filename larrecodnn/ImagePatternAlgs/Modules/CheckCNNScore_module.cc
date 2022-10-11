////////////////////////////////////////////////////////////////////////
// Class:       CheckCNNScore
// Plugin Type: analyzer (art v3_02_06)
// File:        CheckCNNScore_module.cc
//
// Generated at Tue Oct  1 11:26:39 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardataobj/RecoBase/Hit.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"

#include <array>
#include <string>
#include <vector>

namespace pdsp {
  class CheckCNNScore;
}

class pdsp::CheckCNNScore : public art::EDAnalyzer {
public:
  explicit CheckCNNScore(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckCNNScore(CheckCNNScore const&) = delete;
  CheckCNNScore(CheckCNNScore&&) = delete;
  CheckCNNScore& operator=(CheckCNNScore const&) = delete;
  CheckCNNScore& operator=(CheckCNNScore&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  // Declare member data here.
  // Input parameters
  art::InputTag fNNetModuleLabel;      // label of the module used for CNN tagging
  art::InputTag fHitsModuleLabel;      // label of hit finder module
  std::vector<std::string> fNNOutputs; // label of network outputs

  TTree* ftree;
  int run;
  int subrun;
  int event;
  std::vector<short> channel;
  std::vector<short> tpc;
  std::vector<short> plane;
  std::vector<short> wire;
  std::vector<double> charge;
  std::vector<double> peakt;
  std::array<std::vector<double>, 4> scores;
};

pdsp::CheckCNNScore::CheckCNNScore(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fNNetModuleLabel(p.get<art::InputTag>("NNetModuleLabel"))
  , fHitsModuleLabel(p.get<art::InputTag>("HitsModuleLabel"))
  , fNNOutputs(p.get<std::vector<std::string>>("NNOutputs"))
{}

void pdsp::CheckCNNScore::analyze(art::Event const& e)
{

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();
  channel.clear();
  tpc.clear();
  plane.clear();
  wire.clear();
  charge.clear();
  peakt.clear();
  scores = {};

  anab::MVAReader<recob::Hit, 4> hitResults(e, fNNetModuleLabel);

  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (e.getByLabel(fHitsModuleLabel, hitListHandle)) art::fill_ptr_vector(hitlist, hitListHandle);

  // loop over hits
  //  for (size_t h = 0; h < hitResults.size(); ++h) {
  for (auto& hit : hitlist) {

    // Get cnn output for hit h
    std::array<float, 4> cnn_out = hitResults.getOutput(hit);

    if (hit->WireID().Plane == 2) {
      channel.push_back(hit->Channel());
      tpc.push_back(hit->WireID().TPC);
      plane.push_back(hit->WireID().Plane);
      wire.push_back(hit->WireID().Wire);
      charge.push_back(hit->Integral());
      peakt.push_back(hit->PeakTime());
      for (size_t i = 0; i < fNNOutputs.size(); ++i) {
        scores[i].push_back(cnn_out[hitResults.getIndex(fNNOutputs[i])]);
      }
      //      std::cout<<hit->WireID().TPC<<" "
      //               <<hit->WireID().Wire<<" "
      //               <<hit->PeakTime()<<" "
      //               <<cnn_out[hitResults.getIndex("el")]<<" "
      //               <<cnn_out[hitResults.getIndex("inel")]<<" "
      //               <<cnn_out[hitResults.getIndex("none")]<<std::endl;
    }
  }
  if (!channel.empty()) ftree->Fill();
}

void pdsp::CheckCNNScore::beginJob()
{
  art::ServiceHandle<art::TFileService> fileServiceHandle;
  ftree = fileServiceHandle->make<TTree>("ftree", "hit info");
  ftree->Branch("run", &run, "run/I");
  ftree->Branch("event", &event, "event/I");
  ftree->Branch("channel", &channel);
  ftree->Branch("tpc", &tpc);
  ftree->Branch("plane", &plane);
  ftree->Branch("wire", &wire);
  ftree->Branch("charge", &charge);
  ftree->Branch("peakt", &peakt);
  for (size_t i = 0; i < size(scores); ++i) {
    ftree->Branch(Form("score_%ld", i), &scores[i]);
  }
}

DEFINE_ART_MODULE(pdsp::CheckCNNScore)
