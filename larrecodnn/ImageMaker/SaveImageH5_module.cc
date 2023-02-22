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
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SaveImageH5(SaveImageH5 const&) = delete;
  SaveImageH5(SaveImageH5&&) = delete;
  SaveImageH5& operator=(SaveImageH5 const&) = delete;
  SaveImageH5& operator=(SaveImageH5&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual ~SaveImageH5() noexcept;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  hep_hpc::hdf5::File hdffile;
  std::function <decltype(dnn::saveImage)> saveImage_; 
  std::string fHDF5FileName;
};


dnn::SaveImageH5::SaveImageH5(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  saveImage_{art::make_tool<decltype(dnn::saveImage)>(p.get<fhicl::ParameterSet>("imageMaker"), "saveImage")},
  fHDF5FileName(p.get<std::string>("HDF5NAME"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

dnn::SaveImageH5::~SaveImageH5() noexcept
{
}

void dnn::SaveImageH5::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  saveImage_(e, hdffile);
}

void dnn::SaveImageH5::beginJob()
{
  // Implementation of optional member function here.
  hdffile = hep_hpc::hdf5::File(fHDF5FileName, H5F_ACC_TRUNC);
}

DEFINE_ART_MODULE(dnn::SaveImageH5)
