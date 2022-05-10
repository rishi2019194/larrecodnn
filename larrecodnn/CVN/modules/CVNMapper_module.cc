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
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"

#include "dunereco/CVN/art/PixelMapProducer.h"
#include "dunereco/CVN/func/PixelMap.h"
#include "dunereco/CVN/func/TrainingData.h"




namespace cvn {

  template <class T, class U> class CVNMapper : public art::EDProducer {
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

    /// Width of pixel map in tdcs
    unsigned short fTdcWidth;

    /// Length of pixel map in wires
    unsigned short fWireLength;

    /// Length of pixel map in wires
    double fTimeResolution;

    /// Special handling for multiple drift volumes
    bool fMultipleDrifts;
    
    /// ADC threshold for calculating charge from wires directly
    double fThreshold;

    /// PixelMapProducer does the work for us
    PixelMapProducer<T, U> fProducer;

  };

  typedef CVNMapper<recob::Hit, cvn::HitHelper> CVNHitMapper;
  typedef CVNMapper<recob::Wire, cvn::WireHelper> CVNWireMapper;
  typedef CVNMapper<sim::SimChannel, cvn::SimChannelHelper> CVNSimMapper;


  //.......................................................................
  template<class T, class U> CVNMapper<T, U>::CVNMapper(fhicl::ParameterSet const& pset): EDProducer{pset},
  fHitsModuleLabel  (pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>    ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short> ("MinClusterHits")),
  fTdcWidth     (pset.get<unsigned short> ("TdcWidth")),
  fWireLength   (pset.get<unsigned short> ("WireLength")),
  fTimeResolution   (pset.get<unsigned short> ("TimeResolution")),
  fThreshold        (pset.get<double>("Threshold")),
  fProducer      (fWireLength, fTdcWidth, fTimeResolution, fThreshold)
  {

    produces< std::vector<cvn::PixelMap>   >(fClusterPMLabel);

  }

  //......................................................................
  template<class T, class U> CVNMapper<T, U>::~CVNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  template<class T, class U> void CVNMapper<T, U>::beginJob()
  {  }

  //......................................................................
  template <class T, class U> void CVNMapper<T, U>::endJob()
  {
  }

  //......................................................................
  template <class T, class U> void CVNMapper<T, U>::produce(art::Event& evt)
  {
    fProducer.SetMultipleDrifts(fMultipleDrifts);

    std::vector< art::Ptr< T > > hitlist;
    auto hitListHandle = evt.getHandle< std::vector< T > >(fHitsModuleLabel);
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



DEFINE_ART_MODULE(cvn::CVNHitMapper)
DEFINE_ART_MODULE(cvn::CVNWireMapper)
DEFINE_ART_MODULE(cvn::CVNSimMapper)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
