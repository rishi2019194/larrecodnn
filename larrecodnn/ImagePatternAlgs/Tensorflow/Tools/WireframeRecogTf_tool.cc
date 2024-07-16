#include "art/Utilities/ToolMacros.h"
#include "larrecodnn/ImagePatternAlgs/Tensorflow/TF/tf_graph.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IWireframeRecog.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "tensorflow/core/public/session.h"

#include <sys/stat.h>

namespace wframerec_tool {

  class WireframeRecogTf : public IWireframeRecog {
  public:
    explicit WireframeRecogTf(const fhicl::ParameterSet& pset);

    std::vector<std::vector<float>> predictWireframeType(
      const std::vector<std::vector<std::vector<short>>>&) const override;

  private:
    std::unique_ptr<tf::Graph> g; // network graph
    std::string fNNetModelFilePath;
    std::vector<std::string> fNNetOutputPattern;
    bool fUseBundle;
  };

  // ------------------------------------------------------
  WireframeRecogTf::WireframeRecogTf(const fhicl::ParameterSet& pset)
  {
    fNNetModelFilePath = pset.get<std::string>("NNetModelFile", "mymodel.pb");
    fUseBundle = pset.get<bool>("UseSavedModelBundle", false);
    fNNetOutputPattern =
      pset.get<std::vector<std::string>>("NNetOutputPattern", {"cnn_output", "dense_3"});
    if ((fNNetModelFilePath.length() > 3) &&
        (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 3, 3, ".pb") == 0) &&
        !fUseBundle) {
      g = tf::Graph::create(
        findFile(fNNetModelFilePath.c_str()).c_str(), fNNetOutputPattern, fUseBundle);
      if (!g) { throw art::Exception(art::errors::Unknown) << "TF model failed."; }
      mf::LogInfo("WireframeRecogTf") << "TF model loaded.";
    }
    else if ((fNNetModelFilePath.length() > 3) && fUseBundle) {
      g = tf::Graph::create(
        findFile(fNNetModelFilePath.c_str()).c_str(), fNNetOutputPattern, fUseBundle);
      if (!g) { throw art::Exception(art::errors::Unknown) << "TF model failed."; }
      mf::LogInfo("WireframeRecogTf") << "TF model loaded.";
    }
    else {
      mf::LogError("WireframeRecogTf") << "File name extension not supported.";
    }

    setupWframeRecRoiParams(pset);
  }

  // ------------------------------------------------------
  std::vector<std::vector<float>> WireframeRecogTf::predictWireframeType(
    const std::vector<std::vector<std::vector<short>>>& wireframes) const
  {
    if (wireframes.empty() || wireframes.front().empty() || wireframes.front().front().empty()) {
      return std::vector<std::vector<float>>();
      //return std::vector<std::vector<float>>(samples,std::vector<float>(2,0.));
    }

    long long int samples = wireframes.size(), rows = wireframes.front().size(),
                  cols = wireframes.front().front().size();
    //std::cout << " !!!! predictWireframeType: samples = " << samples << ", rows = " << rows << ", cols = " << cols << std::endl;

    std::vector<tensorflow::Tensor> _x;
    _x.push_back(
      tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({samples, rows, cols, 1})));
    auto input_map = _x[0].tensor<float, 4>();
    for (long long int s = 0; s < samples; ++s) {
      const auto& wframe = wireframes[s];
      for (long long int r = 0; r < rows; ++r) {
        const auto& row = wframe[r];
        for (long long int c = 0; c < cols; ++c) {
          input_map(s, r, c, 0) = float(row[c]);
        }
      }
    }
    return g->runx(_x);
  }

}
DEFINE_ART_CLASS_TOOL(wframerec_tool::WireframeRecogTf)
