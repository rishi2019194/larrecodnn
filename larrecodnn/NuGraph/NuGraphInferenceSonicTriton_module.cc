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
#include <chrono>
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

// Function to print elements of a vector<float>
void printVector(const std::vector<float>& vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i];
      // Print space unless it's the last element
      if (i != vec.size() - 1) {
          std::cout << " ";
      }
  }
  std::cout << std::endl;
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
  std::string inference_url;
  bool inference_ssl;
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
  , inference_url(p.get<std::string>("url"))
  , inference_ssl(p.get<bool>("ssl"))
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


  //Here the input should be sent to Triton
  std::string fTritonModelName = "nugraph2";
  std::string fTritonURL = inference_url;
  bool fTritonVerbose = false;
  bool fTritonSSL = inference_ssl;
  std::string fTritonModelVersion = "";
  unsigned fTritonTimeout = 0;
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

  auto hit_table_hit_id_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();
  auto hit_table_local_plane_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();
  auto hit_table_local_time_ptr = std::make_shared<lartriton::TritonInput<float>>();
  auto hit_table_local_wire_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();
  auto hit_table_integral_ptr = std::make_shared<lartriton::TritonInput<float>>();
  auto hit_table_rms_ptr = std::make_shared<lartriton::TritonInput<float>>();
  auto spacepoint_table_spacepoint_id_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();
  auto spacepoint_table_hit_id_u_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();
  auto spacepoint_table_hit_id_v_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();
  auto spacepoint_table_hit_id_y_ptr = std::make_shared<lartriton::TritonInput<int32_t>>();

  hit_table_hit_id_ptr->reserve(1);
  hit_table_local_plane_ptr->reserve(1);
  hit_table_local_time_ptr->reserve(1);
  hit_table_local_wire_ptr->reserve(1);
  hit_table_integral_ptr->reserve(1);
  hit_table_rms_ptr->reserve(1);
  spacepoint_table_spacepoint_id_ptr->reserve(1);
  spacepoint_table_hit_id_u_ptr->reserve(1);
  spacepoint_table_hit_id_v_ptr->reserve(1);
  spacepoint_table_hit_id_y_ptr->reserve(1);

  auto& hit_table_hit_id = hit_table_hit_id_ptr->emplace_back(); 
  auto& hit_table_local_plane = hit_table_local_plane_ptr->emplace_back();
  auto& hit_table_local_time = hit_table_local_time_ptr->emplace_back();
  auto& hit_table_local_wire = hit_table_local_wire_ptr->emplace_back();
  auto& hit_table_integral = hit_table_integral_ptr->emplace_back();
  auto& hit_table_rms = hit_table_rms_ptr->emplace_back();
  auto& spacepoint_table_spacepoint_id = spacepoint_table_spacepoint_id_ptr->emplace_back();
  auto& spacepoint_table_hit_id_u = spacepoint_table_hit_id_u_ptr->emplace_back();
  auto& spacepoint_table_hit_id_v = spacepoint_table_hit_id_v_ptr->emplace_back();
  auto& spacepoint_table_hit_id_y = spacepoint_table_hit_id_y_ptr->emplace_back();

  auto& inputs = triton_client->input();
  for (auto& input_pair : inputs) {
      const std::string& key = input_pair.first;
      auto& triton_input = input_pair.second;   

      if(key == "hit_table_hit_id"){
        for (size_t i = 0; i < hit_table_hit_id_data.size(); ++i) {
            hit_table_hit_id.push_back(hit_table_hit_id_data[i]); 
        }
        triton_input.toServer(hit_table_hit_id_ptr);
      }
      else if(key == "hit_table_local_plane"){
        for (size_t i = 0; i < hit_table_local_plane_data.size(); ++i) {
            hit_table_local_plane.push_back(hit_table_local_plane_data[i]); 
        }
        triton_input.toServer(hit_table_local_plane_ptr);
      }
      else if(key == "hit_table_local_time"){
        for (size_t i = 0; i < hit_table_local_time_data.size(); ++i) {
            hit_table_local_time.push_back(hit_table_local_time_data[i]); 
        }
        triton_input.toServer(hit_table_local_time_ptr);
      }
      else if(key == "hit_table_local_wire"){
        for (size_t i = 0; i < hit_table_local_wire_data.size(); ++i) {
            hit_table_local_wire.push_back(hit_table_local_wire_data[i]); 
        }
        triton_input.toServer(hit_table_local_wire_ptr);
      }
      else if(key == "hit_table_integral"){
        for (size_t i = 0; i < hit_table_integral_data.size(); ++i) {
            hit_table_integral.push_back(hit_table_integral_data[i]); 
        }
        triton_input.toServer(hit_table_integral_ptr);
      }
      else if(key == "hit_table_rms"){
        for (size_t i = 0; i < hit_table_rms_data.size(); ++i) {
            hit_table_rms.push_back(hit_table_rms_data[i]); 
        }
        triton_input.toServer(hit_table_rms_ptr);
      }
      else if(key == "spacepoint_table_spacepoint_id"){
        for (size_t i = 0; i < spacepoint_table_spacepoint_id_data.size(); ++i) {
            spacepoint_table_spacepoint_id.push_back(spacepoint_table_spacepoint_id_data[i]); 
        }
        triton_input.toServer(spacepoint_table_spacepoint_id_ptr);
      }
      else if(key == "spacepoint_table_hit_id_u"){
        for (size_t i = 0; i < spacepoint_table_hit_id_u_data.size(); ++i) {
            spacepoint_table_hit_id_u.push_back(spacepoint_table_hit_id_u_data[i]); 
        }
        triton_input.toServer(spacepoint_table_hit_id_u_ptr);
      }
      else if(key == "spacepoint_table_hit_id_v"){
        for (size_t i = 0; i < spacepoint_table_hit_id_v_data.size(); ++i) {
            spacepoint_table_hit_id_v.push_back(spacepoint_table_hit_id_v_data[i]); 
        }
        triton_input.toServer(spacepoint_table_hit_id_v_ptr);
      }
      else if(key == "spacepoint_table_hit_id_y"){
        for (size_t i = 0; i < spacepoint_table_hit_id_y_data.size(); ++i) {
            spacepoint_table_hit_id_y.push_back(spacepoint_table_hit_id_y_data[i]); 
        }
        triton_input.toServer(spacepoint_table_hit_id_y_ptr);
      }
  }

  auto start = std::chrono::high_resolution_clock::now();
  // ~~~~ Send inference request
  triton_client->dispatch();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time taken for inference: " << elapsed.count() << " seconds" << std::endl;


  // ~~~~ Retrieve inference results
  const auto& triton_output0 = triton_client->output().at("x_semantic_u");
  const auto& prob0 = triton_output0.fromServer<float>();
  size_t triton_input0_elements = std::distance(prob0[0].begin(), prob0[0].end());

  const auto& triton_output1 = triton_client->output().at("x_semantic_v");
  const auto& prob1 = triton_output1.fromServer<float>();
  size_t triton_input1_elements = std::distance(prob1[0].begin(), prob1[0].end());

  const auto& triton_output2 = triton_client->output().at("x_semantic_y");
  const auto& prob2 = triton_output2.fromServer<float>();
  size_t triton_input2_elements = std::distance(prob2[0].begin(), prob2[0].end());

  const auto& triton_output3 = triton_client->output().at("x_filter_u");
  const auto& prob3 = triton_output3.fromServer<float>();
  size_t triton_input3_elements = std::distance(prob3[0].begin(), prob3[0].end());

  const auto& triton_output4 = triton_client->output().at("x_filter_v");
  const auto& prob4 = triton_output4.fromServer<float>();
  size_t triton_input4_elements = std::distance(prob4[0].begin(), prob4[0].end());

  const auto& triton_output5 = triton_client->output().at("x_filter_y");
  const auto& prob5 = triton_output5.fromServer<float>();
  size_t triton_input5_elements = std::distance(prob5[0].begin(), prob5[0].end());
  
  // putting in the resp output vectors
  std::vector<float> x_semantic_u_data;
  x_semantic_u_data.reserve(triton_input0_elements);
  x_semantic_u_data.insert(x_semantic_u_data.end(), prob0[0].begin(), prob0[0].end());

  std::vector<float> x_semantic_v_data;
  x_semantic_v_data.reserve(triton_input1_elements);
  x_semantic_v_data.insert(x_semantic_v_data.end(), prob1[0].begin(), prob1[0].end());

  std::vector<float> x_semantic_y_data;
  x_semantic_y_data.reserve(triton_input2_elements);
  x_semantic_y_data.insert(x_semantic_y_data.end(), prob2[0].begin(), prob2[0].end());

  std::vector<float> x_filter_u_data;
  x_filter_u_data.reserve(triton_input3_elements);
  x_filter_u_data.insert(x_filter_u_data.end(), prob3[0].begin(), prob3[0].end());

  std::vector<float> x_filter_v_data;
  x_filter_v_data.reserve(triton_input4_elements);
  x_filter_v_data.insert(x_filter_v_data.end(), prob4[0].begin(), prob4[0].end());

  std::vector<float> x_filter_y_data;
  x_filter_y_data.reserve(triton_input5_elements);
  x_filter_y_data.insert(x_filter_y_data.end(), prob5[0].begin(), prob5[0].end());

  std::cout<<"Triton Input: "<<std::endl;

  std::cout<<"x_semantic_u: "<<std::endl;
  printVector(x_semantic_u_data);

  std::cout<<"x_semantic_v: "<<std::endl;
  printVector(x_semantic_v_data);

  std::cout<<"x_semantic_y: "<<std::endl;
  printVector(x_semantic_y_data);

  std::cout<<"x_filter_u: "<<std::endl;
  printVector(x_filter_u_data);

  std::cout<<"x_filter_v: "<<std::endl;
  printVector(x_filter_v_data);

  std::cout<<"x_filter_y: "<<std::endl;
  printVector(x_filter_y_data);

  // writing the outputs to the output root file
  if (semanticDecoder) {
    size_t n_cols = 5;
    for (size_t p = 0; p < planes.size(); p++) {
      torch::Tensor s;
      torch::TensorOptions options = torch::TensorOptions().dtype(torch::kFloat32);
      if(planes[p]=="u"){
        size_t n_rows = x_semantic_u_data.size() / n_cols;
        s = torch::from_blob(x_semantic_u_data.data(), {static_cast<int64_t>(n_rows), static_cast<int64_t>(n_cols)}, options);
      }
      else if(planes[p]=="v"){
        size_t n_rows = x_semantic_v_data.size() / n_cols;
        s = torch::from_blob(x_semantic_v_data.data(), {static_cast<int64_t>(n_rows), static_cast<int64_t>(n_cols)}, options);
      }
      else if(planes[p]=="y"){
        size_t n_rows = x_semantic_y_data.size() / n_cols;
        s = torch::from_blob(x_semantic_y_data.data(), {static_cast<int64_t>(n_rows), static_cast<int64_t>(n_cols)}, options);
      }
      else{
        std::cout<<"Error!!"<<std::endl;
      }

      for (int i = 0; i < s.sizes()[0]; ++i) {
        size_t idx = idsmap[p][i];
        std::array<float, 5> input({s[i][0].item<float>(),
                                    s[i][1].item<float>(),
                                    s[i][2].item<float>(),
                                    s[i][3].item<float>(),
                                    s[i][4].item<float>()});
        softmax(input);
        FeatureVector<5> semt = FeatureVector<5>(input);
        (*semtcol)[idx] = semt;
      }
    }
  }
  if (filterDecoder) {
    for (size_t p = 0; p < planes.size(); p++) {
      torch::Tensor f;
      torch::TensorOptions options = torch::TensorOptions().dtype(torch::kFloat32);
      if(planes[p]=="u"){
        int64_t num_elements = x_filter_u_data.size();
        f = torch::from_blob(x_filter_u_data.data(), {num_elements}, options);
      }
      else if(planes[p]=="v"){
        int64_t num_elements = x_filter_v_data.size();
        f = torch::from_blob(x_filter_v_data.data(), {num_elements}, options);
      }
      else if(planes[p]=="y"){
        int64_t num_elements = x_filter_y_data.size();
        f = torch::from_blob(x_filter_y_data.data(), {num_elements}, options);
      }
      else{
        std::cout<<"error!"<<std::endl;
      }

      for (int i = 0; i < f.numel(); ++i) {
        size_t idx = idsmap[p][i];
        std::array<float, 1> input({f[i].item<float>()});
        (*filtcol)[idx] = FeatureVector<1>(input);
      }
    }
  }

  if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
  if (semanticDecoder) {
    e.put(std::move(semtcol), "semantic");
    e.put(std::move(semtdes), "semantic");
  }
  if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
}

DEFINE_ART_MODULE(NuGraphInferenceSonicTriton)
