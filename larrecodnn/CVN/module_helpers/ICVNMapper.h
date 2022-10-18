////////////////////////////////////////////////////////////////////////
// \file    ICVNMapper_module.cc
// \brief   Producer module for creating CVN PixelMap objects
// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef CVN_ICVNMAPPER_H
#define CVN_ICVNMAPPER_H



// C/C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"
#include "larrecodnn/CVN/func/PixelMap.h"

namespace cvn {

  template <class T, class U> class ICVNMapper : public art::EDProducer {
  public:
    explicit ICVNMapper(fhicl::ParameterSet const& pset);
    ~ICVNMapper();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();

  protected:
    /// Module lablel for input clusters
    std::string    fHitsModuleLabel;

    /// Instance lablel for cluster pixelmaps
    std::string    fClusterPMLabel;

    /// Minimum number of hits for cluster to be converted to pixel map
    unsigned short fMinClusterHits;

    /// PixelMapProducer does the work for us
    T fProducer;

  };

}
#endif
