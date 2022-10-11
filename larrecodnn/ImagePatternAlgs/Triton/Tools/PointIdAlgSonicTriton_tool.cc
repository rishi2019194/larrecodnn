////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlgSonicTriton_tool
// Authors:     M.Wang
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/TritonClient.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IPointIdAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <string>
#include <vector>

namespace PointIdAlgTools {

  class PointIdAlgSonicTriton : public IPointIdAlg {
  public:
    explicit PointIdAlgSonicTriton(fhicl::Table<Config> const& table);

    std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) const override;
    std::vector<std::vector<float>> Run(std::vector<std::vector<std::vector<float>>> const& inps,
                                        int samples = -1) const override;

  private:
    std::string fTritonModelName;
    std::string fTritonURL;
    bool fTritonVerbose;
    std::string fTritonModelVersion;
    unsigned fTritonTimeout;
    unsigned fTritonAllowedTries;

    std::unique_ptr<lartriton::TritonClient> triton_client;
  };

  // ------------------------------------------------------
  PointIdAlgSonicTriton::PointIdAlgSonicTriton(fhicl::Table<Config> const& table)
    : img::DataProviderAlg(table())
  {
    // ... Get common config vars
    fNNetOutputs = table().NNetOutputs();
    fPatchSizeW = table().PatchSizeW();
    fPatchSizeD = table().PatchSizeD();
    fCurrentWireIdx = 99999;
    fCurrentScaledDrift = 99999;

    // ... Get "optional" config vars specific to tRTis interface
    fTritonModelName = table().TritonModelName();
    fTritonURL = table().TritonURL();
    fTritonVerbose = table().TritonVerbose();
    fTritonModelVersion = table().TritonModelVersion();
    fTritonAllowedTries = table().TritonAllowedTries();

    // ... Create parameter set for Triton inference client
    fhicl::ParameterSet TritonPset;
    TritonPset.put("serverURL", fTritonURL);
    TritonPset.put("verbose", fTritonVerbose);
    TritonPset.put("modelName", fTritonModelName);
    TritonPset.put("modelVersion", fTritonModelVersion);
    TritonPset.put("timeout", fTritonTimeout);
    TritonPset.put("allowedTries", fTritonAllowedTries);
    TritonPset.put("outputs", "[]");

    // ... Create the Triton inference client
    triton_client = std::make_unique<lartriton::TritonClient>(TritonPset);

    mf::LogInfo("PointIdAlgSonicTriton") << "url: " << fTritonURL;
    mf::LogInfo("PointIdAlgSonicTriton") << "model name: " << fTritonModelName;
    mf::LogInfo("PointIdAlgSonicTriton") << "model version: " << fTritonModelVersion;
    mf::LogInfo("PointIdAlgSonicTriton") << "verbose: " << fTritonVerbose;

    mf::LogInfo("PointIdAlgSonicTriton") << "tensorRT inference context created.";

    resizePatch();
  }

  // ------------------------------------------------------
  std::vector<float> PointIdAlgSonicTriton::Run(std::vector<std::vector<float>> const& inp2d) const
  {
    size_t nrows = inp2d.size();

    triton_client->setBatchSize(1); // set batch size

    // ~~~~ Initialize the inputs
    auto& triton_input = triton_client->input().begin()->second;

    auto data1 = std::make_shared<lartriton::TritonInput<float>>();
    data1->reserve(1);

    // ~~~~ Prepare image for sending to server
    auto& img = data1->emplace_back();
    // ..first flatten the 2d array into contiguous 1d block
    for (size_t ir = 0; ir < nrows; ++ir) {
      img.insert(img.end(), inp2d[ir].begin(), inp2d[ir].end());
    }

    triton_input.toServer(data1); // convert to server format

    // ~~~~ Send inference request
    triton_client->dispatch();

    // ~~~~ Retrieve inference results
    const auto& triton_output0 = triton_client->output().at("em_trk_none_netout/Softmax");
    const auto& prob0 = triton_output0.fromServer<float>();
    auto ncat0 = triton_output0.sizeDims();

    const auto& triton_output1 = triton_client->output().at("michel_netout/Sigmoid");
    const auto& prob1 = triton_output1.fromServer<float>();
    auto ncat1 = triton_output1.sizeDims();

    std::vector<float> out;
    out.reserve(ncat0 + ncat1);
    out.insert(out.end(), prob0[0].begin(), prob0[0].end());
    out.insert(out.end(), prob1[0].begin(), prob1[0].end());

    triton_client->reset();

    return out;
  }

  // ------------------------------------------------------
  std::vector<std::vector<float>> PointIdAlgSonicTriton::Run(
    std::vector<std::vector<std::vector<float>>> const& inps,
    int samples) const
  {
    if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty()) {
      return std::vector<std::vector<float>>();
    }

    if ((samples == -1) || (samples > (long long int)inps.size())) { samples = inps.size(); }

    size_t usamples = samples;
    size_t nrows = inps.front().size();

    triton_client->setBatchSize(usamples); // set batch size

    // ~~~~ Initialize the inputs
    auto& triton_input = triton_client->input().begin()->second;

    auto data1 = std::make_shared<lartriton::TritonInput<float>>();
    data1->reserve(usamples);

    // ~~~~ For each sample, prepare images for sending to server
    for (size_t idx = 0; idx < usamples; ++idx) {
      auto& img = data1->emplace_back();
      // ..first flatten the 2d array into contiguous 1d block
      for (size_t ir = 0; ir < nrows; ++ir) {
        img.insert(img.end(), inps[idx][ir].begin(), inps[idx][ir].end());
      }
    }
    triton_input.toServer(data1); // convert to server format

    // ~~~~ Send inference request
    triton_client->dispatch();

    // ~~~~ Retrieve inference results
    const auto& triton_output0 = triton_client->output().at("em_trk_none_netout/Softmax");
    const auto& prob0 = triton_output0.fromServer<float>();
    auto ncat0 = triton_output0.sizeDims();

    const auto& triton_output1 = triton_client->output().at("michel_netout/Sigmoid");
    const auto& prob1 = triton_output1.fromServer<float>();
    auto ncat1 = triton_output1.sizeDims();

    std::vector<std::vector<float>> out;
    out.reserve(usamples);
    for (unsigned i = 0; i < usamples; i++) {
      out.emplace_back();
      auto& img = out.back();
      img.reserve(ncat0 + ncat1);
      img.insert(img.end(), prob0[i].begin(), prob0[i].end());
      img.insert(img.end(), prob1[i].begin(), prob1[i].end());
    }

    triton_client->reset();

    return out;
  }

}
DEFINE_ART_CLASS_TOOL(PointIdAlgTools::PointIdAlgSonicTriton)
