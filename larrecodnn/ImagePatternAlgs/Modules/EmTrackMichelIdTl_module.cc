/////////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackMichelIdTl
// Module Type: producer
// File:        EmTrackMichelIdTl_module.cc
// Authors:     D.Stefan (dorota.stefan@cern.ch), from DUNE, CERN/NCBJ
//              P.Plonski (pplonski86@gmail.com), from DUNE, WUT
//              R.Sulej (robert.sulej@cern.ch),   from DUNE, FNAL/NCBJ
//              M.Wang (mwang@fnal.gov),          from DUNE, FNAL, tool interface ver
// Module applies CNN to 2D image made of deconvoluted wire waveforms in order
// to distinguish EM-like activity from track-like objects. In addition the activity
// of Michel electrons is recognized. New clusters of hits are produced to include
// also unclustered hits and tag all the activity as EM/track in the same way.
// Note: Michel electrons are best tagged on the level of single hits, since
// clustering may have introduced mistakes (hits of muon and electron clustered
// together).
//
/////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "larrecodnn/ImagePatternAlgs/Modules/EmTrack.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace nnet {

  class EmTrackMichelIdTl : public art::EDProducer {
  public:
    using Parameters = art::EDProducer::Table<EmTrack<4>::Config>;
    explicit EmTrackMichelIdTl(Parameters const& p);

    EmTrackMichelIdTl(EmTrackMichelIdTl const&) = delete;
    EmTrackMichelIdTl(EmTrackMichelIdTl&&) = delete;
    EmTrackMichelIdTl& operator=(EmTrackMichelIdTl const&) = delete;
    EmTrackMichelIdTl& operator=(EmTrackMichelIdTl&&) = delete;

  private:
    void produce(art::Event& e) override;
    EmTrack<4> fEmTrack;
  };
  // ------------------------------------------------------

  EmTrackMichelIdTl::EmTrackMichelIdTl(EmTrackMichelIdTl::Parameters const& p)
    : EDProducer{p}
    , fEmTrack{p(), p.get_PSet().get<std::string>("module_label"), producesCollector()}
  {}
  // ------------------------------------------------------

  void EmTrackMichelIdTl::produce(art::Event& evt)
  {
    fEmTrack.produce(evt);
  }
  // ------------------------------------------------------

  DEFINE_ART_MODULE(EmTrackMichelIdTl)

}
