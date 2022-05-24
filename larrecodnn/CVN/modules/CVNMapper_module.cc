////////////////////////////////////////////////////////////////////////
// \file    CVNMapper_module.cc
// \brief   Producer module for creating CVN PixelMap objects
// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

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

  class CVNMapper : public art::EDProducer {
  public:
    explicit CVNMapper(fhicl::ParameterSet const& pset);
    ~CVNMapper();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();

  private:
    /// Module lablel for input clusters
    std::string    fHitsModuleLabel;

    /// Instance lablel for cluster pixelmaps
    std::string    fClusterPMLabel;

    /// Minimum number of hits for cluster to be converted to pixel map
    unsigned short fMinClusterHits;

    /// PixelMapProducer does the work for us
    PixelMapHitProducer fProducer;

  };

  //.......................................................................
  CVNMapper::CVNMapper(fhicl::ParameterSet const& pset): EDProducer{pset},
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>    ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short> ("MinClusterHits")),
  fProducer(pset.get<fhicl::ParameterSet> ("PixelMapProducer"))
  {

    produces< std::vector<cvn::PixelMap> >(fClusterPMLabel);

  }

  //......................................................................
  CVNMapper::~CVNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNMapper::beginJob()
  {  }

  //......................................................................
  void CVNMapper::endJob()
  {
  }

  //......................................................................
  void CVNMapper::produce(art::Event& evt)
  {

    std::vector< art::Ptr< recob::Hit > > hitlist;
    auto hitListHandle = evt.getHandle< std::vector< recob::Hit > >(fHitsModuleLabel);
    if (hitListHandle)
      art::fill_ptr_vector(hitlist, hitListHandle);

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::PixelMap> >
      pmCol(new std::vector<cvn::PixelMap>);

    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    PixelMap pm = fProducer.CreateMap(detProp, hitlist);
    auto nhits = fProducer.NROI();
    pm.SetTotHits(nhits);
   
    if (nhits > fMinClusterHits) 
      pmCol->push_back(pm);
    
    evt.put(std::move(pmCol), fClusterPMLabel);
  }

  //----------------------------------------------------------------------

DEFINE_ART_MODULE(cvn::CVNMapper)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
