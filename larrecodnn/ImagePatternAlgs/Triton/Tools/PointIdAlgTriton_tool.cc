////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlgTriton_tool
// Authors:     M.Wang,                                         FNAL, 2021: Nvidia Triton inf client
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IPointIdAlg.h"

// Nvidia Triton inference server client includes
#include "grpc_client.h"

namespace ni = nvidia::inferenceserver;
namespace nic = nvidia::inferenceserver::client;

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace PointIdAlgTools {

  class PointIdAlgTriton : public IPointIdAlg {
  public:
    explicit PointIdAlgTriton(fhicl::Table<Config> const& table);

    std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) const override;
    std::vector<std::vector<float>> Run(std::vector<std::vector<std::vector<float>>> const& inps,
                                        int samples = -1) const override;

  private:
    std::string fTritonModelName;
    std::string fTritonURL;
    bool fTritonVerbose;
    std::string fTritonModelVersion;

    std::unique_ptr<nic::InferenceServerGrpcClient> triton_client;
    inference::ModelMetadataResponse triton_modmet;
    inference::ModelConfigResponse triton_modcfg;
    mutable std::vector<int64_t> triton_inpshape;
    nic::InferOptions triton_options;
  };

  // ------------------------------------------------------
  PointIdAlgTriton::PointIdAlgTriton(fhicl::Table<Config> const& table)
    : img::DataProviderAlg(table()), triton_options("")
  {
    // ... Get common config vars
    fNNetOutputs = table().NNetOutputs();
    fPatchSizeW = table().PatchSizeW();
    fPatchSizeD = table().PatchSizeD();
    fCurrentWireIdx = 99999;
    fCurrentScaledDrift = 99999;

    // ... Get "optional" config vars specific to Triton interface
    fTritonModelName = table().TritonModelName();
    fTritonURL = table().TritonURL();
    fTritonVerbose = table().TritonVerbose();
    fTritonModelVersion = table().TritonModelVersion();

    // ... Create the Triton inference client
    auto err = nic::InferenceServerGrpcClient::Create(&triton_client, fTritonURL, fTritonVerbose);
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "error: unable to create client for inference: " << err << std::endl;
    }

    // ... Get the model metadata and config information
    err = triton_client->ModelMetadata(&triton_modmet, fTritonModelName, fTritonModelVersion);
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "error: failed to get model metadata: " << err << std::endl;
    }
    err = triton_client->ModelConfig(&triton_modcfg, fTritonModelName, fTritonModelVersion);
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "error: failed to get model config: " << err << std::endl;
    }

    // ... Set up shape vector needed when creating inference input
    triton_inpshape.push_back(1); // initialize batch_size to 1
    triton_inpshape.push_back(triton_modmet.inputs(0).shape(1));
    triton_inpshape.push_back(triton_modmet.inputs(0).shape(2));
    triton_inpshape.push_back(triton_modmet.inputs(0).shape(3));

    // ... Set up Triton inference client options
    triton_options.model_name_ = fTritonModelName;
    triton_options.model_version_ = fTritonModelVersion;

    mf::LogInfo("PointIdAlgTriton") << "url: " << fTritonURL;
    mf::LogInfo("PointIdAlgTriton") << "model name: " << fTritonModelName;
    mf::LogInfo("PointIdAlgTriton") << "model version: " << fTritonModelVersion;
    mf::LogInfo("PointIdAlgTriton") << "verbose: " << fTritonVerbose;

    mf::LogInfo("PointIdAlgTriton") << "tensorRT inference context created.";

    resizePatch();
  }

  // ------------------------------------------------------
  std::vector<float> PointIdAlgTriton::Run(std::vector<std::vector<float>> const& inp2d) const
  {
    size_t nrows = inp2d.size(), ncols = inp2d.front().size();

    triton_inpshape.at(0) = 1; // set batch size

    // ~~~~ Initialize the inputs

    nic::InferInput* triton_input;
    auto err = nic::InferInput::Create(&triton_input,
                                       triton_modmet.inputs(0).name(),
                                       triton_inpshape,
                                       triton_modmet.inputs(0).datatype());
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton") << "unable to get input: " << err << std::endl;
    }
    std::shared_ptr<nic::InferInput> triton_input_ptr(triton_input);
    std::vector<nic::InferInput*> triton_inputs = {triton_input_ptr.get()};

    // ~~~~ Register the mem address of 1st byte of image and #bytes in image

    err = triton_input_ptr->Reset();
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "failed resetting Triton model input: " << err << std::endl;
    }

    size_t sbuff_byte_size = (nrows * ncols) * sizeof(float);
    std::vector<float> fa(sbuff_byte_size);

    // ..flatten the 2d array into contiguous 1d block
    for (size_t ir = 0; ir < nrows; ++ir) {
      std::copy(inp2d[ir].begin(), inp2d[ir].end(), fa.begin() + (ir * ncols));
    }
    err = triton_input_ptr->AppendRaw(reinterpret_cast<uint8_t*>(fa.data()), sbuff_byte_size);
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "failed setting Triton input: " << err << std::endl;
    }

    // ~~~~ Send inference request

    nic::InferResult* results;

    err = triton_client->Infer(&results, triton_options, triton_inputs);
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "failed sending Triton synchronous infer request: " << err << std::endl;
    }
    std::shared_ptr<nic::InferResult> results_ptr;
    results_ptr.reset(results);

    // ~~~~ Retrieve inference results

    std::vector<float> out;

    const float* prb0;
    size_t rbuff0_byte_size; // size of result buffer in bytes
    results_ptr->RawData(
      triton_modmet.outputs(0).name(), (const uint8_t**)&prb0, &rbuff0_byte_size);
    size_t ncat0 = rbuff0_byte_size / sizeof(float);

    const float* prb1;
    size_t rbuff1_byte_size; // size of result buffer in bytes
    results_ptr->RawData(
      triton_modmet.outputs(1).name(), (const uint8_t**)&prb1, &rbuff1_byte_size);
    size_t ncat1 = rbuff1_byte_size / sizeof(float);

    for (unsigned j = 0; j < ncat0; j++)
      out.push_back(*(prb0 + j));
    for (unsigned j = 0; j < ncat1; j++)
      out.push_back(*(prb1 + j));

    return out;
  }

  // ------------------------------------------------------
  std::vector<std::vector<float>> PointIdAlgTriton::Run(
    std::vector<std::vector<std::vector<float>>> const& inps,
    int samples) const
  {
    if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty()) {
      return std::vector<std::vector<float>>();
    }

    if ((samples == -1) || (samples > (long long int)inps.size())) { samples = inps.size(); }

    size_t usamples = samples;
    size_t nrows = inps.front().size(), ncols = inps.front().front().size();

    triton_inpshape.at(0) = usamples; // set batch size

    // ~~~~ Initialize the inputs

    nic::InferInput* triton_input;
    auto err = nic::InferInput::Create(&triton_input,
                                       triton_modmet.inputs(0).name(),
                                       triton_inpshape,
                                       triton_modmet.inputs(0).datatype());
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton") << "unable to get input: " << err << std::endl;
    }
    std::shared_ptr<nic::InferInput> triton_input_ptr(triton_input);
    std::vector<nic::InferInput*> triton_inputs = {triton_input_ptr.get()};

    // ~~~~ For each sample, register the mem address of 1st byte of image and #bytes in image
    err = triton_input_ptr->Reset();
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "failed resetting Triton model input: " << err << std::endl;
    }

    size_t sbuff_byte_size = (nrows * ncols) * sizeof(float);
    std::vector<std::vector<float>> fa(usamples, std::vector<float>(sbuff_byte_size));

    for (size_t idx = 0; idx < usamples; ++idx) {
      // ..first flatten the 2d array into contiguous 1d block
      for (size_t ir = 0; ir < nrows; ++ir) {
        std::copy(inps[idx][ir].begin(), inps[idx][ir].end(), fa[idx].begin() + (ir * ncols));
      }
      err =
        triton_input_ptr->AppendRaw(reinterpret_cast<uint8_t*>(fa[idx].data()), sbuff_byte_size);
      if (!err.IsOk()) {
        throw cet::exception("PointIdAlgTriton")
          << "failed setting Triton input: " << err << std::endl;
      }
    }

    // ~~~~ Send inference request

    nic::InferResult* results;

    err = triton_client->Infer(&results, triton_options, triton_inputs);
    if (!err.IsOk()) {
      throw cet::exception("PointIdAlgTriton")
        << "failed sending Triton synchronous infer request: " << err << std::endl;
    }
    std::shared_ptr<nic::InferResult> results_ptr;
    results_ptr.reset(results);

    // ~~~~ Retrieve inference results

    std::vector<std::vector<float>> out;

    const float* prb0;
    size_t rbuff0_byte_size; // size of result buffer in bytes
    results_ptr->RawData(
      triton_modmet.outputs(0).name(), (const uint8_t**)&prb0, &rbuff0_byte_size);
    size_t ncat0 = rbuff0_byte_size / (usamples * sizeof(float));

    const float* prb1;
    size_t rbuff1_byte_size; // size of result buffer in bytes
    results_ptr->RawData(
      triton_modmet.outputs(1).name(), (const uint8_t**)&prb1, &rbuff1_byte_size);
    size_t ncat1 = rbuff1_byte_size / (usamples * sizeof(float));

    for (unsigned i = 0; i < usamples; i++) {
      std::vector<float> vprb;
      for (unsigned j = 0; j < ncat0; j++)
        vprb.push_back(*(prb0 + i * ncat0 + j));
      for (unsigned j = 0; j < ncat1; j++)
        vprb.push_back(*(prb1 + i * ncat1 + j));
      out.push_back(vprb);
    }

    return out;
  }

}
DEFINE_ART_CLASS_TOOL(PointIdAlgTools::PointIdAlgTriton)
