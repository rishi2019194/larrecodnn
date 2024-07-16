#include "art/Utilities/ToolMacros.h"
#include "larrecodnn/ImagePatternAlgs/Tensorflow/TF/tf_graph.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IWaveformDenoise.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "tensorflow/core/public/session.h"

#include <sys/stat.h>

namespace wavdenoise_tool {

  class WaveformDenoiseTf : public IWaveformDenoise {
  public:
    explicit WaveformDenoiseTf(const fhicl::ParameterSet& pset);

    std::vector<std::vector<float>> applyDenoisingAE(
      const std::vector<std::vector<float>>&) const override;

  private:
    std::unique_ptr<tf::Graph> g; // network graph
    std::string fNNetModelFilePath;
    std::vector<std::string> fNNetOutputPattern;
    bool fUseBundle;
  };

  // ------------------------------------------------------
  WaveformDenoiseTf::WaveformDenoiseTf(const fhicl::ParameterSet& pset)
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
      mf::LogInfo("WaveformDenoiseTf") << "TF model loaded (pb format).";
    }
    else if ((fNNetModelFilePath.length() > 3) && fUseBundle) {
      g = tf::Graph::create(
        findFile(fNNetModelFilePath.c_str()).c_str(), fNNetOutputPattern, fUseBundle);
      if (!g) { throw art::Exception(art::errors::Unknown) << "TF model failed."; }
      mf::LogInfo("WaveformDenoiseTf") << "TF model loaded (SavedModel format).";
    }
    else {
      mf::LogError("WaveformDenoiseTf") << "File name extension not supported.";
    }

    setupWaveDenoiseParams(pset);
  }

  // ------------------------------------------------------
  std::vector<std::vector<float>> WaveformDenoiseTf::applyDenoisingAE(
    const std::vector<std::vector<float>>& waveforms) const
  {
    if (waveforms.empty() || waveforms.front().empty()) {
      return std::vector<std::vector<float>>();
    }

    long long int samples = waveforms.size(), numtcks = waveforms.front().size();

    //std::cout<<"Samples: "<<samples<<", Ticks: "<<numtcks<<std::endl;
    std::vector<tensorflow::Tensor> _x;
    _x.push_back(
      tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({samples, numtcks, 1})));
    auto input_map = _x[0].tensor<float, 3>();
    for (long long int s = 0; s < samples; ++s) {
      const auto& wvfrm = waveforms[s];
      for (long long int t = 0; t < numtcks; ++t) {
        input_map(s, t, 0) = wvfrm[t];
      }
    }

    return g->runae(_x);
  }

}
DEFINE_ART_CLASS_TOOL(wavdenoise_tool::WaveformDenoiseTf)
