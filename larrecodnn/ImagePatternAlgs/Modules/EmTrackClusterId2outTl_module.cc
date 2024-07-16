/////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackClusterId2outTl
// Module Type: producer
// File:        EmTrackClusterId2outTl_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// Module applies CNN to 2D image made of deconvoluted wire waveforms in order
// to distinguish EM-like activity from track-like objects. New clusters of
// hits are produced to include also unclustered hits and tag everything in
// a common way.
// NOTE: This module uses 2-output CNN models, see EmTrackClusterId and
// EmTrackMichelClusterId for usage of 3 and 4-output models.
//
/////////////////////////////////////////////////////////////////////////////////

#include "larrecodnn/ImagePatternAlgs/Modules/EmTrack.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace nnet {

  class EmTrackClusterId2outTl : public art::EDProducer {
  public:
    using Parameters = art::EDProducer::Table<EmTrack<2>::Config>;
    explicit EmTrackClusterId2outTl(Parameters const& p);

    EmTrackClusterId2outTl(EmTrackClusterId2outTl const&) = delete;
    EmTrackClusterId2outTl(EmTrackClusterId2outTl&&) = delete;
    EmTrackClusterId2outTl& operator=(EmTrackClusterId2outTl const&) = delete;
    EmTrackClusterId2outTl& operator=(EmTrackClusterId2outTl&&) = delete;

  private:
    void produce(art::Event& e) override;
    EmTrack<2> fEmTrack;
  };
  // ------------------------------------------------------

  EmTrackClusterId2outTl::EmTrackClusterId2outTl(EmTrackClusterId2outTl::Parameters const& p)
    : EDProducer{p}
    , fEmTrack{p(), p.get_PSet().get<std::string>("module_label"), producesCollector()}
  {}
  // ------------------------------------------------------

  void EmTrackClusterId2outTl::produce(art::Event& evt)
  {
    fEmTrack.produce(evt);
  }
  // ------------------------------------------------------

  DEFINE_ART_MODULE(EmTrackClusterId2outTl)

}
