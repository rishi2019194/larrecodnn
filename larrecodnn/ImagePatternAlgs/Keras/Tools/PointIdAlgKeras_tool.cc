////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlgKeras_tool (tool version of Keras model interface in PointIdAlg)
// Authors:     D.Stefan (Dorota.Stefan@ncbj.gov.pl),      from DUNE,   CERN/NCBJ, since May 2016
//              R.Sulej (Robert.Sulej@cern.ch),            from DUNE,   FNAL/NCBJ, since May 2016
//              P.Plonski,                                 from DUNE,   WUT,       since May 2016
//              D.Smith,                                   from LArIAT, BU, 2017: real data dump
//              M.Wang,                                    from DUNE,   FNAL, 2020: tool version
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"

#include "larrecodnn/ImagePatternAlgs/Keras/keras_model.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IPointIdAlg.h"

#include <sys/stat.h>

namespace PointIdAlgTools {

  class PointIdAlgKeras : public IPointIdAlg {
  public:
    explicit PointIdAlgKeras(const fhicl::ParameterSet& pset)
      : PointIdAlgKeras(fhicl::Table<Config>(pset, {})())
    {}
    explicit PointIdAlgKeras(const Config& config);

    std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) const override;
    std::vector<std::vector<float>> Run(std::vector<std::vector<std::vector<float>>> const& inps,
                                        int samples = -1) const override;

  private:
    std::unique_ptr<keras::KerasModel> m;
    std::string fNNetModelFilePath;
    std::string findFile(const char* fileName) const;
  };

  // ------------------------------------------------------
  PointIdAlgKeras::PointIdAlgKeras(const Config& config) : img::DataProviderAlg(config)
  {
    // ... Get common config vars
    fNNetOutputs = config.NNetOutputs();
    fPatchSizeW = config.PatchSizeW();
    fPatchSizeD = config.PatchSizeD();
    fCurrentWireIdx = 99999;
    fCurrentScaledDrift = 99999;

    // ... Get "optional" config vars specific to tf interface
    std::string s_cfgvr;
    if (config.NNetModelFile(s_cfgvr)) { fNNetModelFilePath = s_cfgvr; }
    else {
      fNNetModelFilePath = "mycnn";
    }

    if ((fNNetModelFilePath.length() > 5) &&
        (fNNetModelFilePath.compare(fNNetModelFilePath.length() - 5, 5, ".nnet") == 0)) {
      m = std::make_unique<keras::KerasModel>(findFile(fNNetModelFilePath.c_str()).c_str());
      mf::LogInfo("PointIdAlgKeras") << "Keras model loaded.";
    }
    else {
      mf::LogError("PointIdAlgKeras") << "File name extension not supported.";
    }

    resizePatch();
  }

  // ------------------------------------------------------
  std::string PointIdAlgKeras::findFile(const char* fileName) const
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
  std::vector<float> PointIdAlgKeras::Run(std::vector<std::vector<float>> const& inp2d) const
  {
    std::vector<std::vector<std::vector<float>>> inp3d;
    inp3d.push_back(inp2d); // lots of copy, should add 2D to keras...

    keras::DataChunk2D sample;
    sample.set_data(inp3d);
    return m->compute_output(&sample);
  }

  // ------------------------------------------------------
  std::vector<std::vector<float>> PointIdAlgKeras::Run(
    std::vector<std::vector<std::vector<float>>> const& inps,
    int samples) const
  {

    if ((samples == 0) || inps.empty() || inps.front().empty() || inps.front().front().empty()) {
      return std::vector<std::vector<float>>();
    }

    if ((samples == -1) || (samples > (long long int)inps.size())) { samples = inps.size(); }

    std::vector<std::vector<float>> out;

    for (long long int s = 0; s < samples; ++s) {
      std::vector<std::vector<std::vector<float>>> inp3d;
      inp3d.push_back(inps[s]); // lots of copy, should add 2D to keras...

      keras::DataChunk* sample = new keras::DataChunk2D();
      sample->set_data(inp3d); // and more copy...
      out.push_back(m->compute_output(sample));
      delete sample;
    }

    return out;
  }

}
DEFINE_ART_CLASS_TOOL(PointIdAlgTools::PointIdAlgKeras)
