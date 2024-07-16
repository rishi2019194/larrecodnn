////////////////////////////////////////////////////////////////////////
// Class:       NuGraphAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuGraphAnalyzer_module.cc
//
// Generated at Mon Nov 20 13:42:17 2023 by Giuseppe Cerati using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// saving output
#include "TTree.h"
#include "art_root_io/TFileService.h"

#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"

class NuGraphAnalyzer;

using std::vector;

class NuGraphAnalyzer : public art::EDAnalyzer {
public:
  explicit NuGraphAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuGraphAnalyzer(NuGraphAnalyzer const&) = delete;
  NuGraphAnalyzer(NuGraphAnalyzer&&) = delete;
  NuGraphAnalyzer& operator=(NuGraphAnalyzer const&) = delete;
  NuGraphAnalyzer& operator=(NuGraphAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // Declare member data here.
  TTree *_treeHit, *_treeEvt;
  int _run, _subrun, _event, _id, _wire, _plane;
  float _x_filter, _MIP, _HIP, _shower, _michel, _diffuse, _time;
  float _vtx_x, _vtx_y, _vtx_z;
};

NuGraphAnalyzer::NuGraphAnalyzer(fhicl::ParameterSet const& p) : EDAnalyzer{p} // ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  _treeHit = tfs->make<TTree>("NuGraphHitOutput", "NuGraphHitOutput");
  _treeHit->Branch("run", &_run, "run/I");
  _treeHit->Branch("subrun", &_subrun, "subrun/I");
  _treeHit->Branch("event", &_event, "event/I");
  _treeHit->Branch("id", &_id, "id/I");
  _treeHit->Branch("wire", &_wire, "wire/I");
  _treeHit->Branch("plane", &_plane, "plane/I");
  _treeHit->Branch("x_filter", &_x_filter, "x_filter/F");
  _treeHit->Branch("MIP", &_MIP, "MIP/F");
  _treeHit->Branch("HIP", &_HIP, "HIP/F");
  _treeHit->Branch("shower", &_shower, "shower/F");
  _treeHit->Branch("michel", &_michel, "michel/F");
  _treeHit->Branch("diffuse", &_diffuse, "diffuse/F");
  _treeHit->Branch("time", &_time, "time/F");
  _treeEvt = tfs->make<TTree>("NuGraphEventOutput", "NuGraphEventOutput");
  _treeEvt->Branch("run", &_run, "run/I");
  _treeEvt->Branch("subrun", &_subrun, "subrun/I");
  _treeEvt->Branch("event", &_event, "event/I");
  _treeEvt->Branch("vtx_x", &_vtx_x, "vtx_x/F");
  _treeEvt->Branch("vtx_y", &_vtx_y, "vtx_y/F");
  _treeEvt->Branch("vtx_z", &_vtx_z, "vtx_z/F");
}

void NuGraphAnalyzer::analyze(art::Event const& e)
{

  auto GNNDescription = e.getHandle<anab::MVADescription<5>>(art::InputTag("NuGraph", "semantic"));

  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
    e,
    GNNDescription->dataTag(), //tag of the hit collection we ran the GNN on
    proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
    proxy::withParallelData<anab::FeatureVector<5>>(art::InputTag("NuGraph", "semantic")));

  std::cout << hitsWithScores.size() << std::endl;
  for (auto& h : hitsWithScores) {
    const auto& assocFilter = h.get<anab::FeatureVector<1>>();
    const auto& assocSemantic = h.get<anab::FeatureVector<5>>();
    _event = e.event();
    _subrun = e.subRun();
    _run = e.run();
    _id = h.index();
    _x_filter = assocFilter.at(0);
    _MIP = assocSemantic.at(GNNDescription->getIndex("MIP"));
    _HIP = assocSemantic.at(GNNDescription->getIndex("HIP"));
    _shower = assocSemantic.at(GNNDescription->getIndex("shower"));
    _michel = assocSemantic.at(GNNDescription->getIndex("michel"));
    _diffuse = assocSemantic.at(GNNDescription->getIndex("diffuse"));
    _treeHit->Fill();
  }

  auto PredVertexColl = e.getHandle<std::vector<recob::Vertex>>(art::InputTag("NuGraph", "vertex"));
  if (PredVertexColl.isValid() && PredVertexColl->size() > 0) { //there should be only one
    _vtx_x = PredVertexColl->at(0).position().X();
    _vtx_y = PredVertexColl->at(0).position().Y();
    _vtx_z = PredVertexColl->at(0).position().Z();
    _treeEvt->Fill();
  }
}

DEFINE_ART_MODULE(NuGraphAnalyzer)
