////////////////////////////////////////////////////////////////////////
// Class:       SaveImageH5
// Plugin Type: analyzer (Unknown Unknown)
// File:        SaveImageH5_module.cc
//
// Generated at Thu Feb  9 17:15:36 2023 by Tingjun Yang using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ImageMaker.h"

namespace dnn {
  class SaveImageH5;
}

class dnn::SaveImageH5 : public art::EDAnalyzer {
public:
  explicit SaveImageH5(fhicl::ParameterSet const& p);

  SaveImageH5(SaveImageH5 const&) = delete;
  SaveImageH5(SaveImageH5&&) = delete;
  SaveImageH5& operator=(SaveImageH5 const&) = delete;
  SaveImageH5& operator=(SaveImageH5&&) = delete;

  void analyze(art::Event const& e) override;
  virtual ~SaveImageH5() noexcept;

  void beginJob() override;

private:
  hep_hpc::hdf5::File hdffile;
  std::unique_ptr<dnn::ImageMaker> saveImage_;
  std::string fHDF5FileName;
};

dnn::SaveImageH5::SaveImageH5(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , saveImage_{art::make_tool<dnn::ImageMaker>(p.get<fhicl::ParameterSet>("imageMaker"))}
  , fHDF5FileName(p.get<std::string>("HDF5NAME"))
{}

dnn::SaveImageH5::~SaveImageH5() noexcept {}

void dnn::SaveImageH5::analyze(art::Event const& e)
{
  saveImage_->saveImage(e, hdffile);
}

void dnn::SaveImageH5::beginJob()
{
  hdffile = hep_hpc::hdf5::File(fHDF5FileName, H5F_ACC_TRUNC);
}

DEFINE_ART_MODULE(dnn::SaveImageH5)
