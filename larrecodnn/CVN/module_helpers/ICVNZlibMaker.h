////////////////////////////////////////////////////////////////////////
// \file    ICVNZlibMaker_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
//           - wrote the zlib code used in this module
////////////////////////////////////////////////////////////////////////
#ifndef CVN_ICVNZLIBMAKER_H
#define CVN_ICVNZLIBMAKER_H

// C/C++ includes
#include <cstdlib>
#include <iostream>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "boost/filesystem.hpp"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"

// Data products
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/TrainingData.h"

// Compression
#include "math.h"
#include "zlib.h"

#include "TH1.h"

namespace fs = boost::filesystem;

namespace cvn {

  class ICVNZlibMaker : public art::EDAnalyzer {
  public:
    explicit ICVNZlibMaker(fhicl::ParameterSet const& pset);
    ~ICVNZlibMaker();

    void beginJob() override;
    void analyze(const art::Event& evt) override {}
    void reconfigure(const fhicl::ParameterSet& pset);

  protected:
    std::string fOutputDir;
    std::string fPixelMapInput;
    bool fSetLog;
    std::vector<bool> fReverseViews;
    unsigned int fPlaneLimit;
    unsigned int fTDCLimit;

    std::string out_dir;
    CVNImageUtils fImage;

    template <class T>
    void write_files(TrainingData<T> td, std::string evtid) = delete;
  };
}
#endif
