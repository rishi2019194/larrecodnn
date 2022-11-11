////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlgTf_tool (tool version of TensorFlow model interface in PointIdAlg)
// Authors:     D.Stefan (Dorota.Stefan@ncbj.gov.pl),      from DUNE,   CERN/NCBJ, since May 2016
//              R.Sulej (Robert.Sulej@cern.ch),            from DUNE,   FNAL/NCBJ, since May 2016
//              P.Plonski,                                 from DUNE,   WUT,       since May 2016
//              D.Smith,                                   from LArIAT, BU, 2017: real data dump
//              M.Wang,                                    from DUNE,   FNAL, 2020: tool version
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"

#include "larrecodnn/ImagePatternAlgs/Tensorflow/TF/tf_graph.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IPointIdAlg.h"
#include "tensorflow/core/public/session.h"

#include <sys/stat.h>

namespace PointIdAlgTools {

  class PointIdAlgTf : public IPointIdAlg {
  public:
    explicit PointIdAlgTf(fhicl::Table<Config> const& table);

    std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) const override;
    std::vector<std::vector<float>> Run(std::vector<std::vector<std::vector<float>>> const& inps,
                                        int samples = -1) const override;

  protected:
    std::string findFile(const char* fileName) const;

  private:
    std::unique_ptr<tf::Graph> g; // network graph
    std::vector<std::string> fNNetOutputPattern;
    std::string fNNetModelFilePath;
  };

  // ------------------------------------------------------
  PointIdAlgTf::PointIdAlgTf(fhicl::Table<Config> const& table) : img::DataProviderAlg(table())
  {
    // ... Get common config vars
    fNNetOutputs = table().NNetOutputs();
    fPatchSizeW = table().PatchSizeW();
    fPatchSizeD = table().PatchSizeD();
    fCurrentWireIdx = 99999;
    fCurrentScaledDrift = 99999;

    // ... Get "optional" config vars specific to tf interface
    std::string s_cfgvr;
    if (table().NNetModelFile(s_cfgvr)) { fNNetModelFilePath = s_cfgvr; }
    else {
      fNNetModelFilePath = "mycnn";
    }
    std::vector<std::string> vs_cfgvr;
    if (table().NNetOutputPattern(vs_cfgvr)) { fNNetOutputPattern = vs_cfgvr; }
    else {
      fNNetOutputPattern = {"cnn_output", "_netout"};
    }

    if ((fNNetModelFilePath.length() > 3) &&
        (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 3, 3, ".pb") == 0)) {
      g = tf::Graph::create(findFile(fNNetModelFilePath.c_str()).c_str(), fNNetOutputPattern);
      if (!g) { throw art::Exception(art::errors::Unknown) << "TF model failed."; }
      mf::LogInfo("PointIdAlgTf") << "TF model loaded.";
    }
    else {
      mf::LogError("PointIdAlgTf") << "File name extension not supported.";
    }

    resizePatch();
  }

  // ------------------------------------------------------
  std::string PointIdAlgTf::findFile(const char* fileName) const
  {
    std::string fname_out;
    cet::search_path sp("FW_SEARCH_PATH");
    if (!sp.find_file(fileName, fname_out)) {
      struct stat buffer;
      if (stat(fileName, &buffer) == 0) { fname_out = fileName; }
      else {
        throw art::Exception(art::errors::NotFound) << "Could not find the model file " << fileName;
      }
    }
    return fname_out;
  }

  // ------------------------------------------------------
  std::vector<float> PointIdAlgTf::Run(std::vector<std::vector<float>> const& inp2d) const
  {
    long long int rows = inp2d.size(), cols = inp2d.front().size();

    std::vector<tensorflow::Tensor> _x;
    _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT,
                                    tensorflow::TensorShape({1, rows, cols, 1})));
    auto input_map = _x[0].tensor<float, 4>();
    for (long long int r = 0; r < rows; ++r) {
      const auto& row = inp2d[r];
      for (long long int c = 0; c < cols; ++c) {
        input_map(0, r, c, 0) = row[c];
      }
    }

    auto out = g->runx(_x);
    if (!out.empty())
      return out.front();
    else
      return std::vector<float>();
  }

  // ------------------------------------------------------
  std::vector<std::vector<float>> PointIdAlgTf::Run(
    std::vector<std::vector<std::vector<float>>> const& inps,
    int samples) const
  {

    if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty()) {
      return std::vector<std::vector<float>>();
    }

    if ((samples == -1) || (samples > (long long int)inps.size())) { samples = inps.size(); }

    long long int rows = inps.front().size(), cols = inps.front().front().size();

    std::vector<tensorflow::Tensor> _x;
    _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT,
                                    tensorflow::TensorShape({samples , rows, cols, 1})));
    auto input_map = _x[0].tensor<float, 4>();
    for (long long int s = 0; s < samples; ++s) {
      const auto& sample = inps[s];
      for (long long int r = 0; r < rows; ++r) {
        const auto& row = sample[r];
        for (long long int c = 0; c < cols; ++c) {
          input_map(s, r, c, 0) = row[c];
        }
      }
    }
    return g->runx(_x);
  }

}
DEFINE_ART_CLASS_TOOL(PointIdAlgTools::PointIdAlgTf)
