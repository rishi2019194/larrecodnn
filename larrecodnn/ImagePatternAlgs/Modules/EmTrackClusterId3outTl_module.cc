/////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackClusterId3outTl
// Module Type: producer
// File:        EmTrackClusterId3outTl_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// Module applies CNN to 2D image made of deconvoluted wire waveforms in order
// to distinguish EM-like activity from track-like objects. New clusters of
// hits are produced to include also unclustered hits and tag everything in
// a common way.
// NOTE: This module uses 3-output CNN models, see EmTrackMichelClusterId for
// usage of 4-output models and EmTrackClusterId3outTl2out_module.cc for 2-output
// models.
//
/////////////////////////////////////////////////////////////////////////////////

#include "larrecodnn/ImagePatternAlgs/Modules/EmTrack.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace nnet {

  class EmTrackClusterId3outTl : public art::EDProducer {
  public:
    using Parameters = art::EDProducer::Table<EmTrack<3>::Config>;
    explicit EmTrackClusterId3outTl(Parameters const& p);

    EmTrackClusterId3outTl(EmTrackClusterId3outTl const&) = delete;
    EmTrackClusterId3outTl(EmTrackClusterId3outTl&&) = delete;
    EmTrackClusterId3outTl& operator=(EmTrackClusterId3outTl const&) = delete;
    EmTrackClusterId3outTl& operator=(EmTrackClusterId3outTl&&) = delete;

  private:
    void produce(art::Event& e) override;
    EmTrack<3> fEmTrack;
  };
  // ------------------------------------------------------

  EmTrackClusterId3outTl::EmTrackClusterId3outTl(EmTrackClusterId3outTl::Parameters const& p)
    : EDProducer{p}
    , fEmTrack{p(), p.get_PSet().get<std::string>("module_label"), producesCollector()}
  {}
  // ------------------------------------------------------

  void EmTrackClusterId3outTl::produce(art::Event& evt)
  {
    fEmTrack.produce(evt);
  }
  // ------------------------------------------------------

  DEFINE_ART_MODULE(EmTrackClusterId3outTl)

}
