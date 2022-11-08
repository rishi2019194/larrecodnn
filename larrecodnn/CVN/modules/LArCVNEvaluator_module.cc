////////////////////////////////////////////////////////////////////////
// \file    LArCVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
//          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/Result.h"
#include "larrecodnn/CVN/interfaces/ITFNetHandler.h"

namespace lcvn {

  class LArCVNEvaluator : public art::EDProducer {
  public:
    explicit LArCVNEvaluator(fhicl::ParameterSet const& pset);
    ~LArCVNEvaluator();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();

  private:
    /// Module label for input pixel maps
    std::string fPixelMapInput;
    std::string fResultLabel;

    /// Can use Caffe or Tensorflow
    std::string fCVNType;

    //lcvn::CaffeNetHandler fCaffeHandler;
    std::unique_ptr<lcvn::ITFNetHandler> fTFHandler;

    /// Number of outputs fron neural net
    //unsigned int fNOutput;

    /// If there are multiple pixel maps per event can we use them?
    bool fMultiplePMs;
  };

  //.......................................................................
  LArCVNEvaluator::LArCVNEvaluator(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fPixelMapInput(pset.get<std::string>("PixelMapInput"))
    , fResultLabel(pset.get<std::string>("ResultLabel"))
    , fCVNType(pset.get<std::string>("CVNType"))
    ,
    //fCaffeHandler       (pset.get<fhicl::ParameterSet> ("CaffeNetHandler")),
    fTFHandler{art::make_tool<ITFNetHandler>(pset.get<fhicl::ParameterSet>("TFHandler"))}
    ,
    //fNOutput       (fCaffeHandler.NOutput()),
    fMultiplePMs(pset.get<bool>("MultiplePMs"))
  {
    produces<std::vector<lcvn::Result>>(fResultLabel);
  }
  //......................................................................
  LArCVNEvaluator::~LArCVNEvaluator()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void LArCVNEvaluator::beginJob() {}

  //......................................................................
  void LArCVNEvaluator::endJob() {}

  //......................................................................
  void LArCVNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr<std::vector<Result>> resultCol(new std::vector<Result>);

    /// Load in the pixel maps
    std::vector<art::Ptr<lcvn::PixelMap>> pixelmaplist;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto pixelmapListHandle = evt.getHandle<std::vector<lcvn::PixelMap>>(itag1);
    if (pixelmapListHandle) art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    if (fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow") {

      // If we have a pixel map then use the TF interface to give us a prediction
      if (pixelmaplist.size() > 0) {

        std::vector<std::vector<float>> networkOutput = fTFHandler->Predict(*pixelmaplist[0]);
        // lcvn::Result can now take a vector of floats and works out the number of outputs
        resultCol->emplace_back(networkOutput);

        // Classify other pixel maps if they exist
        if (fMultiplePMs) {
          for (unsigned int p = 1; p < pixelmaplist.size(); ++p) {
            std::vector<std::vector<float>> output = fTFHandler->Predict(*pixelmaplist[p]);
            resultCol->emplace_back(output);
          }
        }
      }
    }
    else {
      mf::LogError("LArCVNEvaluator::produce")
        << "CVN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("LArCVNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }

    evt.put(std::move(resultCol), fResultLabel);
  }

  DEFINE_ART_MODULE(lcvn::LArCVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
