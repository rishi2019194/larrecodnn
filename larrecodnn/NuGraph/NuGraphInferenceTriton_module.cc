////////////////////////////////////////////////////////////////////////
// Class:       NuGraphInferenceTriton
// Plugin Type: producer (Unknown Unknown)
// File:        NuGraphInferenceTriton_module.cc
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

#include <array>
#include <limits>
#include <memory>

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

class NuGraphInferenceTriton;

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


class NuGraphInferenceTriton : public art::EDProducer {
public:
  explicit NuGraphInferenceTriton(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  NuGraphInferenceTriton(NuGraphInferenceTriton const&) = delete;
  NuGraphInferenceTriton(NuGraphInferenceTriton&&) = delete;
  NuGraphInferenceTriton& operator=(NuGraphInferenceTriton const&) = delete;
  NuGraphInferenceTriton& operator=(NuGraphInferenceTriton&&) = delete;

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

NuGraphInferenceTriton::NuGraphInferenceTriton(fhicl::ParameterSet const& p)
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

void NuGraphInferenceTriton::produce(art::Event& e)
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
  bool verbose = false;
  std::string url(inference_url);
  tc::Headers http_headers;
  uint32_t client_timeout = 0;
  bool use_ssl = inference_ssl;
  std::string root_certificates;
  std::string private_key;
  std::string certificate_chain;
  grpc_compression_algorithm compression_algorithm =
      grpc_compression_algorithm::GRPC_COMPRESS_NONE;
  bool test_use_cached_channel = false;
  bool use_cached_channel = true;

  // the element-wise difference.
  std::string model_name = "nugraph2";
  std::string model_version = "";

  // Create a InferenceServerGrpcClient instance to communicate with the
  // server using gRPC protocol.
  std::unique_ptr<tc::InferenceServerGrpcClient> client;
  tc::SslOptions ssl_options = tc::SslOptions();
  std::string err;
  if (use_ssl) {
    ssl_options.root_certificates = root_certificates;
    ssl_options.private_key = private_key;
    ssl_options.certificate_chain = certificate_chain;
    err = "unable to create secure grpc client";
  } else {
    err = "unable to create grpc client";
  }
  // Run with the same name to ensure cached channel is not used
  int numRuns = test_use_cached_channel ? 2 : 1;
  for (int i = 0; i < numRuns; ++i) {
    FAIL_IF_ERR(
        tc::InferenceServerGrpcClient::Create(
            &client, url, verbose, use_ssl, ssl_options, tc::KeepAliveOptions(),
            use_cached_channel),
        err);



    std::vector<int64_t> hit_table_shape{int64_t(hit_table_hit_id_data.size())};
    std::vector<int64_t> spacepoint_table_shape{int64_t(spacepoint_table_spacepoint_id_data.size())};
  
    // Initialize the inputs with the data.
    tc::InferInput* hit_table_hit_id;
    tc::InferInput* hit_table_local_plane;
    tc::InferInput* hit_table_local_time;
    tc::InferInput* hit_table_local_wire;
    tc::InferInput* hit_table_integral;
    tc::InferInput* hit_table_rms;

    tc::InferInput* spacepoint_table_spacepoint_id;
    tc::InferInput* spacepoint_table_hit_id_u;
    tc::InferInput* spacepoint_table_hit_id_v;
    tc::InferInput* spacepoint_table_hit_id_y;


    FAIL_IF_ERR(
        tc::InferInput::Create(&hit_table_hit_id, "hit_table_hit_id", hit_table_shape, "INT32"),
        "unable to get hit_table_hit_id");
    std::shared_ptr<tc::InferInput> hit_table_hit_id_ptr;
    hit_table_hit_id_ptr.reset(hit_table_hit_id);

    FAIL_IF_ERR(
        tc::InferInput::Create(&hit_table_local_plane, "hit_table_local_plane", hit_table_shape, "INT32"),
        "unable to get hit_table_local_plane");
    std::shared_ptr<tc::InferInput> hit_table_local_plane_ptr;
    hit_table_local_plane_ptr.reset(hit_table_local_plane);

    FAIL_IF_ERR(
        tc::InferInput::Create(&hit_table_local_time, "hit_table_local_time", hit_table_shape, "FP32"),
        "unable to get hit_table_local_time");
    std::shared_ptr<tc::InferInput> hit_table_local_time_ptr;
    hit_table_local_time_ptr.reset(hit_table_local_time);

    FAIL_IF_ERR(
        tc::InferInput::Create(&hit_table_local_wire, "hit_table_local_wire", hit_table_shape, "INT32"),
        "unable to get hit_table_local_wire");
    std::shared_ptr<tc::InferInput> hit_table_local_wire_ptr;
    hit_table_local_wire_ptr.reset(hit_table_local_wire);

    FAIL_IF_ERR(
        tc::InferInput::Create(&hit_table_integral, "hit_table_integral", hit_table_shape, "FP32"),
        "unable to get hit_table_integral");
    std::shared_ptr<tc::InferInput> hit_table_integral_ptr;
    hit_table_integral_ptr.reset(hit_table_integral);

    FAIL_IF_ERR(
        tc::InferInput::Create(&hit_table_rms, "hit_table_rms", hit_table_shape, "FP32"),
        "unable to get hit_table_rms");
    std::shared_ptr<tc::InferInput> hit_table_rms_ptr;
    hit_table_rms_ptr.reset(hit_table_rms);


    FAIL_IF_ERR(
        tc::InferInput::Create(&spacepoint_table_spacepoint_id, "spacepoint_table_spacepoint_id", spacepoint_table_shape, "INT32"),
        "unable to get spacepoint_table_spacepoint_id");
    std::shared_ptr<tc::InferInput> spacepoint_table_spacepoint_id_ptr;
    spacepoint_table_spacepoint_id_ptr.reset(spacepoint_table_spacepoint_id);

    FAIL_IF_ERR(
        tc::InferInput::Create(&spacepoint_table_hit_id_u, "spacepoint_table_hit_id_u", spacepoint_table_shape, "INT32"),
        "unable to get spacepoint_table_spacepoint_hit_id_u");
    std::shared_ptr<tc::InferInput> spacepoint_table_hit_id_u_ptr;
    spacepoint_table_hit_id_u_ptr.reset(spacepoint_table_hit_id_u);

    FAIL_IF_ERR(
        tc::InferInput::Create(&spacepoint_table_hit_id_v, "spacepoint_table_hit_id_v", spacepoint_table_shape, "INT32"),
        "unable to get spacepoint_table_spacepoint_hit_id_v");
    std::shared_ptr<tc::InferInput> spacepoint_table_hit_id_v_ptr;
    spacepoint_table_hit_id_v_ptr.reset(spacepoint_table_hit_id_v);

    FAIL_IF_ERR(
        tc::InferInput::Create(&spacepoint_table_hit_id_y, "spacepoint_table_hit_id_y", spacepoint_table_shape, "INT32"),
        "unable to get spacepoint_table_spacepoint_hit_id_y");
    std::shared_ptr<tc::InferInput> spacepoint_table_hit_id_y_ptr;
    spacepoint_table_hit_id_y_ptr.reset(spacepoint_table_hit_id_y);



    FAIL_IF_ERR(
        hit_table_hit_id_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&hit_table_hit_id_data[0]),
            hit_table_hit_id_data.size() * sizeof(float)),
        "unable to set data for hit_table_hit_id");

    FAIL_IF_ERR(
        hit_table_local_plane_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&hit_table_local_plane_data[0]),
            hit_table_local_plane_data.size() * sizeof(float)),
        "unable to set data for hit_table_local_plane");
    
    FAIL_IF_ERR(
        hit_table_local_time_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&hit_table_local_time_data[0]),
            hit_table_local_time_data.size() * sizeof(float)),
        "unable to set data for hit_table_local_time");

    FAIL_IF_ERR(
        hit_table_local_wire_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&hit_table_local_wire_data[0]),
            hit_table_local_wire_data.size() * sizeof(float)),
        "unable to set data for hit_table_local_wire");

    FAIL_IF_ERR(
        hit_table_integral_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&hit_table_integral_data[0]),
            hit_table_integral_data.size() * sizeof(float)),
        "unable to set data for hit_table_integral");

    FAIL_IF_ERR(
        hit_table_rms_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&hit_table_rms_data[0]),
            hit_table_rms_data.size() * sizeof(float)),
        "unable to set data for hit_table_rms");

    FAIL_IF_ERR(
        spacepoint_table_spacepoint_id_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&spacepoint_table_spacepoint_id_data[0]),
            spacepoint_table_spacepoint_id_data.size() * sizeof(float)),
        "unable to set data for spacepoint_table_spacepoint_id");

    FAIL_IF_ERR(
        spacepoint_table_hit_id_u_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&spacepoint_table_hit_id_u_data[0]),
            spacepoint_table_hit_id_u_data.size() * sizeof(float)),
        "unable to set data for spacepoint_table_hit_id_u");

    FAIL_IF_ERR(
        spacepoint_table_hit_id_v_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&spacepoint_table_hit_id_v_data[0]),
            spacepoint_table_hit_id_v_data.size() * sizeof(float)),
        "unable to set data for spacepoint_table_hit_id_v");

    FAIL_IF_ERR(
        spacepoint_table_hit_id_y_ptr->AppendRaw(
            reinterpret_cast<uint8_t*>(&spacepoint_table_hit_id_y_data[0]),
            spacepoint_table_hit_id_y_data.size() * sizeof(float)),
        "unable to set data for spacepoint_table_hit_id_y");

    // Generate the outputs to be requested.
    tc::InferRequestedOutput* x_semantic_u;
    tc::InferRequestedOutput* x_semantic_v;
    tc::InferRequestedOutput* x_semantic_y;
    tc::InferRequestedOutput* x_filter_u;
    tc::InferRequestedOutput* x_filter_v;
    tc::InferRequestedOutput* x_filter_y;


    FAIL_IF_ERR(
        tc::InferRequestedOutput::Create(&x_semantic_u, "x_semantic_u"),
        "unable to get 'x_semantic_u'");
    std::shared_ptr<tc::InferRequestedOutput> x_semantic_u_ptr;
    x_semantic_u_ptr.reset(x_semantic_u);

    FAIL_IF_ERR(
        tc::InferRequestedOutput::Create(&x_semantic_v, "x_semantic_v"),
        "unable to get 'x_semantic_v'");
    std::shared_ptr<tc::InferRequestedOutput> x_semantic_v_ptr;
    x_semantic_v_ptr.reset(x_semantic_v);

    FAIL_IF_ERR(
        tc::InferRequestedOutput::Create(&x_semantic_y, "x_semantic_y"),
        "unable to get 'x_semantic_y'");
    std::shared_ptr<tc::InferRequestedOutput> x_semantic_y_ptr;
    x_semantic_y_ptr.reset(x_semantic_y);

    FAIL_IF_ERR(
        tc::InferRequestedOutput::Create(&x_filter_u, "x_filter_u"),
        "unable to get 'x_filter_u'");
    std::shared_ptr<tc::InferRequestedOutput> x_filter_u_ptr;
    x_filter_u_ptr.reset(x_filter_u);

    FAIL_IF_ERR(
        tc::InferRequestedOutput::Create(&x_filter_v, "x_filter_v"),
        "unable to get 'x_filter_v'");
    std::shared_ptr<tc::InferRequestedOutput> x_filter_v_ptr;
    x_filter_v_ptr.reset(x_filter_v);

    FAIL_IF_ERR(
        tc::InferRequestedOutput::Create(&x_filter_y, "x_filter_y"),
        "unable to get 'x_filter_y'");
    std::shared_ptr<tc::InferRequestedOutput> x_filter_y_ptr;
    x_filter_y_ptr.reset(x_filter_y);


    // The inference settings. Will be using default for now.
    tc::InferOptions options(model_name);
    options.model_version_ = model_version;
    options.client_timeout_ = client_timeout;

    std::cout<<options.model_name_<<std::endl;
    std::cout<<options.model_version_<<std::endl;
    std::cout<<options.client_timeout_<<std::endl;

    std::vector<tc::InferInput*> inputs = {hit_table_hit_id_ptr.get(), hit_table_local_plane_ptr.get(), hit_table_local_time_ptr.get(), \
                                          hit_table_local_wire_ptr.get(), hit_table_integral_ptr.get(), hit_table_rms_ptr.get(), \
                                          spacepoint_table_spacepoint_id_ptr.get(), spacepoint_table_hit_id_u_ptr.get(), \
                                          spacepoint_table_hit_id_v_ptr.get(), spacepoint_table_hit_id_y_ptr.get()};

    // Iterate over the vector and print each InferInput object's value
    for (size_t i = 0; i < inputs.size(); ++i) {
        // Assuming tc::InferInput has a method or member variable to get its value
        // Replace 'getValue()' with the appropriate method or member access to get the value.
        std::cout << "Element " << inputs[i]->Name() << ": " << inputs[i]->Shape() << std::endl;
    }

    std::vector<const tc::InferRequestedOutput*> outputs = {
        x_semantic_u_ptr.get(), x_semantic_v_ptr.get(), \
        x_semantic_y_ptr.get(), x_filter_u_ptr.get(), x_filter_v_ptr.get(), \
        x_filter_y_ptr.get()};

    // Iterate over the vector and print each InferInput object's value
    for (size_t i = 0; i < outputs.size(); ++i) {
        // Assuming tc::InferInput has a method or member variable to get its value
        // Replace 'getValue()' with the appropriate method or member access to get the value.
        std::cout << "Element " << outputs[i]->Name() << std::endl;
    }

    tc::InferResult* results;
    FAIL_IF_ERR(
        client->Infer(
            &results, options, inputs, outputs, http_headers,
            compression_algorithm),
        "unable to run model");
    std::shared_ptr<tc::InferResult> results_ptr;
    results_ptr.reset(results);

    // Get pointers to the result returned...


    const float* x_semantic_u_data;
    size_t x_semantic_u_byte_size;
    FAIL_IF_ERR(
        results_ptr->RawData(
            "x_semantic_u", (const uint8_t**)&x_semantic_u_data, &x_semantic_u_byte_size),
        "unable to get result data for 'x_semantic_u'");

    const float* x_semantic_v_data;
    size_t x_semantic_v_byte_size;
    FAIL_IF_ERR(
        results_ptr->RawData(
            "x_semantic_v", (const uint8_t**)&x_semantic_v_data, &x_semantic_v_byte_size),
        "unable to get result data for 'x_semantic_v'");

    const float* x_semantic_y_data;
    size_t x_semantic_y_byte_size;
    FAIL_IF_ERR(
        results_ptr->RawData(
            "x_semantic_y", (const uint8_t**)&x_semantic_y_data, &x_semantic_y_byte_size),
        "unable to get result data for 'x_semantic_y'");

    const float* x_filter_u_data;
    size_t x_filter_u_byte_size;
    FAIL_IF_ERR(
        results_ptr->RawData(
            "x_filter_u", (const uint8_t**)&x_filter_u_data, &x_filter_u_byte_size),
        "unable to get result data for 'x_filter_u'");

    const float* x_filter_v_data;
    size_t x_filter_v_byte_size;
    FAIL_IF_ERR(
        results_ptr->RawData(
            "x_filter_v", (const uint8_t**)&x_filter_v_data, &x_filter_v_byte_size),
        "unable to get result data for 'x_filter_v'");

    const float* x_filter_y_data;
    size_t x_filter_y_byte_size;
    FAIL_IF_ERR(
        results_ptr->RawData(
            "x_filter_y", (const uint8_t**)&x_filter_y_data, &x_filter_y_byte_size),
        "unable to get result data for 'x_filter_y'");



    std::cout<<"Trition output: "<<std::endl;  

    std::cout<<"x_semantic_u: "<<std::endl;
    printFloatArray(x_semantic_u_data, x_semantic_u_byte_size/sizeof(float));
    std::cout << std::endl;

    std::cout<<"x_semantic_v: "<<std::endl;
    printFloatArray(x_semantic_v_data, x_semantic_v_byte_size/sizeof(float));
    std::cout << std::endl;

    std::cout<<"x_semantic_y: "<<std::endl;
    printFloatArray(x_semantic_y_data, x_semantic_y_byte_size/sizeof(float));
    std::cout << std::endl;

    std::cout<<"x_filter_u: "<<std::endl;
    printFloatArray(x_filter_u_data, x_filter_u_byte_size/sizeof(float));
    std::cout << std::endl;

    std::cout<<"x_filter_v: "<<std::endl;
    printFloatArray(x_filter_v_data, x_filter_v_byte_size/sizeof(float));
    std::cout << std::endl;

    std::cout<<"x_filter_y: "<<std::endl;
    printFloatArray(x_filter_y_data, x_filter_y_byte_size/sizeof(float));
    std::cout << std::endl;
  
    if (semanticDecoder) {
      size_t n_cols = 5;
      for (size_t p = 0; p < planes.size(); p++) {
        torch::Tensor s;
        torch::TensorOptions options = torch::TensorOptions().dtype(torch::kFloat32);
        if(planes[p]=="u"){
          size_t n_rows = x_semantic_u_byte_size / (n_cols*sizeof(float));
          s = torch::from_blob((void*)x_semantic_u_data, {static_cast<int64_t>(n_rows), static_cast<int64_t>(n_cols)}, options);
        }
        else if(planes[p]=="v"){
          size_t n_rows = x_semantic_v_byte_size / (n_cols*sizeof(float));
          s = torch::from_blob((void*)x_semantic_v_data, {static_cast<int64_t>(n_rows), static_cast<int64_t>(n_cols)}, options);
        }
        else if(planes[p]=="y"){
          size_t n_rows = x_semantic_y_byte_size / (n_cols*sizeof(float));
          s = torch::from_blob((void*)x_semantic_y_data, {static_cast<int64_t>(n_rows), static_cast<int64_t>(n_cols)}, options);
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
          int64_t num_elements = x_filter_u_byte_size/sizeof(float);
          f = torch::from_blob((void*)x_filter_u_data, {num_elements}, options);
        }
        else if(planes[p]=="v"){
          int64_t num_elements = x_filter_v_byte_size/sizeof(float);
          f = torch::from_blob((void*)x_filter_v_data, {num_elements}, options);
        }
        else if(planes[p]=="y"){
          int64_t num_elements = x_filter_y_byte_size/sizeof(float);
          f = torch::from_blob((void*)x_filter_y_data, {num_elements}, options);
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
  }
  if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
  if (semanticDecoder) {
    e.put(std::move(semtcol), "semantic");
    e.put(std::move(semtdes), "semantic");
  }
  if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
}

DEFINE_ART_MODULE(NuGraphInferenceTriton)
