////////////////////////////////////////////////////////////////////////
// Class:       NuGraphInferenceSonicTriton
// Plugin Type: producer (Unknown Unknown)
// File:        NuGraphInferenceSonicTriton_module.cc
//
// Generated at Tue Nov 14 14:41:30 2023 by Giuseppe Cerati using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/types/Table.h"

#include <array>
#include <limits>
#include <memory>

#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/TritonClient.h"
#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/TritonData.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h" //this creates a conflict with torch script if included before it...

#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <torch/torch.h>

#include "grpc_client.h"

class NuGraphInferenceSonicTriton;

using anab::FeatureVector;
using anab::MVADescription;
using recob::Hit;
using recob::SpacePoint;
using std::array;
using std::vector;


#define FAIL_IF_ERR(X, MSG)                                        \
  {                                                                \
    tc::Error err = (X);                                           \
    if (!err.IsOk()) {                                             \
      std::cerr << "error: " << (MSG) << ": " << err << std::endl; \
      exit(1);                                                     \
    }                                                              \
  }
namespace tc = triton::client;

namespace {

  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

  template <typename T, size_t N>
  void softmax(std::array<T, N>& arr)
  {
    T m = -std::numeric_limits<T>::max();
    for (size_t i = 0; i < arr.size(); i++) {
      if (arr[i] > m) { m = arr[i]; }
    }
    T sum = 0.0;
    for (size_t i = 0; i < arr.size(); i++) {
      sum += expf(arr[i] - m);
    }
    T offset = m + logf(sum);
    for (size_t i = 0; i < arr.size(); i++) {
      arr[i] = expf(arr[i] - offset);
    }
    return;
  }
}

// Function to convert string to integer
int stoi(const std::string &str) {
    std::istringstream iss(str);
    int num;
    iss >> num;
    return num;
}

void printFloatArray(const float* data, size_t num_elements) {
    for (size_t i = 0; i < num_elements; ++i) {
        std::cout << data[i];
        if (i < num_elements - 1) {
            std::cout << " ";
        }
    }
    std::cout << std::endl;
}


class NuGraphInferenceSonicTriton : public art::EDProducer {
public:
  explicit NuGraphInferenceSonicTriton(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  NuGraphInferenceSonicTriton(NuGraphInferenceSonicTriton const&) = delete;
  NuGraphInferenceSonicTriton(NuGraphInferenceSonicTriton&&) = delete;
  NuGraphInferenceSonicTriton& operator=(NuGraphInferenceSonicTriton const&) = delete;
  NuGraphInferenceSonicTriton& operator=(NuGraphInferenceSonicTriton&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  vector<std::string> planes;
  art::InputTag hitInput;
  art::InputTag spsInput;
  size_t minHits;
  bool debug;
  // vector<vector<float>> avgs;
  // vector<vector<float>> devs;
  bool filterDecoder;
  bool semanticDecoder;
  bool vertexDecoder;
};

NuGraphInferenceSonicTriton::NuGraphInferenceSonicTriton(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , planes(p.get<vector<std::string>>("planes"))
  , hitInput(p.get<art::InputTag>("hitInput"))
  , spsInput(p.get<art::InputTag>("spsInput"))
  , minHits(p.get<size_t>("minHits"))
  , debug(p.get<bool>("debug"))
  , filterDecoder(p.get<bool>("filterDecoder"))
  , semanticDecoder(p.get<bool>("semanticDecoder"))
  , vertexDecoder(p.get<bool>("vertexDecoder"))
{

  // for (size_t ip = 0; ip < planes.size(); ++ip) {
  //   avgs.push_back(p.get<vector<float>>("avgs_" + planes[ip]));
  //   devs.push_back(p.get<vector<float>>("devs_" + planes[ip]));
  // }

  if (filterDecoder) { produces<vector<FeatureVector<1>>>("filter"); }
  //
  if (semanticDecoder) {
    produces<vector<FeatureVector<5>>>("semantic");
    produces<MVADescription<5>>("semantic");
  }
  //
  if (vertexDecoder) { produces<vector<recob::Vertex>>("vertex"); }
}

void NuGraphInferenceSonicTriton::produce(art::Event& e)
{

  art::Handle<vector<Hit>> hitListHandle;
  vector<art::Ptr<Hit>> hitlist;
  if (e.getByLabel(hitInput, hitListHandle)) { art::fill_ptr_vector(hitlist, hitListHandle); }

  std::unique_ptr<vector<FeatureVector<1>>> filtcol(
    new vector<FeatureVector<1>>(hitlist.size(), FeatureVector<1>(std::array<float, 1>({-1.}))));

  std::unique_ptr<vector<FeatureVector<5>>> semtcol(new vector<FeatureVector<5>>(
    hitlist.size(), FeatureVector<5>(std::array<float, 5>({-1., -1., -1., -1., -1.}))));
  std::unique_ptr<MVADescription<5>> semtdes(
    new MVADescription<5>(hitListHandle.provenance()->moduleLabel(),
                          "semantic",
                          {"MIP", "HIP", "shower", "michel", "diffuse"}));

  std::unique_ptr<vector<recob::Vertex>> vertcol(new vector<recob::Vertex>());

  if (debug) std::cout << "Hits size=" << hitlist.size() << std::endl;
  if (hitlist.size() < minHits) {
    if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
    if (semanticDecoder) {
      e.put(std::move(semtcol), "semantic");
      e.put(std::move(semtdes), "semantic");
    }
    if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
    return;
  }
  vector<vector<size_t>> idsmap(planes.size(), vector<size_t>());
  vector<size_t> idsmapRev(hitlist.size(), hitlist.size());
  for (auto h : hitlist) {
    idsmap[h->View()].push_back(h.key());
    idsmapRev[h.key()] = idsmap[h->View()].size() - 1;
  }

  // event id
  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  array<int, 3> evtID;
  evtID[0] = run;
  evtID[1] = subrun;
  evtID[2] = event;  

  // hit table
  vector<int32_t> hit_table_hit_id_data;
  vector<int32_t> hit_table_local_plane_data;
  vector<float>  hit_table_local_time_data;
  vector<int32_t> hit_table_local_wire_data;
  vector<float>  hit_table_integral_data;
  vector<float>  hit_table_rms_data;
  for (auto h : hitlist) {
    hit_table_hit_id_data.push_back(h.key());
    hit_table_local_plane_data.push_back(h->View());
    hit_table_local_time_data.push_back(h->PeakTime());
    hit_table_local_wire_data.push_back(h->WireID().Wire);
    hit_table_integral_data.push_back(h->Integral());
    hit_table_rms_data.push_back(h->RMS());
  }


  // Get spacepoints from the event record
  art::Handle<vector<SpacePoint>> spListHandle;
  vector<art::Ptr<SpacePoint>> splist;
  if (e.getByLabel(spsInput, spListHandle)) { art::fill_ptr_vector(splist, spListHandle); }
  // Get assocations from spacepoints to hits
  vector<vector<art::Ptr<Hit>>> sp2Hit(splist.size());
  if (splist.size() > 0) {
    art::FindManyP<Hit> fmp(spListHandle, e, "sps");
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    }
  }

  // space point table
  vector<int32_t> spacepoint_table_spacepoint_id_data;
  vector<int32_t> spacepoint_table_hit_id_u_data;
  vector<int32_t> spacepoint_table_hit_id_v_data;
  vector<int32_t> spacepoint_table_hit_id_y_data;
  for (size_t i = 0; i < splist.size(); ++i) {
    spacepoint_table_spacepoint_id_data.push_back(i);
    spacepoint_table_hit_id_u_data.push_back(-1);
    spacepoint_table_hit_id_v_data.push_back(-1);
    spacepoint_table_hit_id_y_data.push_back(-1);
    for (size_t j = 0; j < sp2Hit[i].size(); ++j) {
      if (sp2Hit[i][j]->View()==0) spacepoint_table_hit_id_u_data.back() = sp2Hit[i][j].key();
      if (sp2Hit[i][j]->View()==1) spacepoint_table_hit_id_v_data.back() = sp2Hit[i][j].key();
      if (sp2Hit[i][j]->View()==2) spacepoint_table_hit_id_y_data.back() = sp2Hit[i][j].key();
    }
  }

  std::cout<<hit_table_hit_id_data.size()<<std::endl;
  std::cout<<hit_table_local_plane_data.size()<<std::endl;
  std::cout<<hit_table_local_time_data.size()<<std::endl;
  std::cout<<hit_table_local_wire_data.size()<<std::endl;
  std::cout<<hit_table_integral_data.size()<<std::endl;
  std::cout<<hit_table_rms_data.size()<<std::endl;
  std::cout<<spacepoint_table_spacepoint_id_data.size()<<std::endl;
  std::cout<<spacepoint_table_hit_id_u_data.size()<<std::endl; 
  std::cout<<spacepoint_table_hit_id_v_data.size()<<std::endl;
  std::cout<<spacepoint_table_hit_id_y_data.size()<<std::endl;


  //Here the input should be sent to Triton
  std::string fTritonModelName = "nugraph2";
  std::string fTritonURL = "triton.fnal.gov:443";
  bool fTritonVerbose = true;
  bool fTritonSSL = true;
  std::string fTritonModelVersion = "";
  unsigned fTritonTimeout = 1000;
  unsigned fTritonAllowedTries = 1;
  std::unique_ptr<lartriton::TritonClient> triton_client;

  // ... Create parameter set for Triton inference client
  fhicl::ParameterSet TritonPset;
  TritonPset.put("serverURL", fTritonURL);
  TritonPset.put("verbose", fTritonVerbose);
  TritonPset.put("ssl", fTritonSSL);
  TritonPset.put("modelName", fTritonModelName);
  TritonPset.put("modelVersion", fTritonModelVersion);
  TritonPset.put("timeout", fTritonTimeout);
  TritonPset.put("allowedTries", fTritonAllowedTries);
  TritonPset.put("outputs", "[]");

  // ... Create the Triton inference client
  triton_client = std::make_unique<lartriton::TritonClient>(TritonPset);

  triton_client->setBatchSize(1); // set batch size

  // ~~~~ Initialize the inputs
  auto& triton_input = triton_client->input().begin()->second;

  // Create data1
  auto data1 = std::make_shared<lartriton::TritonInput<float>>();
  data1->reserve(1); 

  // Emplace back hit_table_hit_id_data into data1
  auto& img = data1->emplace_back(); 

  // Fill img with elements from hit_table_hit_id_data
  for (size_t i = 0; i < hit_table_hit_id_data.size(); ++i) {
      img.push_back(hit_table_hit_id_data[i]); 
  }
  std::cout<<data1->size()<<std::endl;
  triton_input.toServer(data1); // convert to server format
  if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
  if (semanticDecoder) {
    e.put(std::move(semtcol), "semantic");
    e.put(std::move(semtdes), "semantic");
  }
  if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
}

DEFINE_ART_MODULE(NuGraphInferenceSonicTriton)
