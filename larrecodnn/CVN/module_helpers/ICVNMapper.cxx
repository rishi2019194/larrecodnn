#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"

namespace cvn {

  template <class T, class U>
  ICVNMapper<T, U>::ICVNMapper(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fHitsModuleLabel(pset.get<std::string>("HitsModuleLabel"))
    , fClusterPMLabel(pset.get<std::string>("ClusterPMLabel"))
    , fMinClusterHits(pset.get<unsigned short>("MinClusterHits"))
    , fProducer(pset.get<fhicl::ParameterSet>("PixelMapProducer"))
  {

    produces<std::vector<cvn::PixelMap>>(fClusterPMLabel);
  }

  //......................................................................
  template <class T, class U>
  ICVNMapper<T, U>::~ICVNMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  template <class T, class U>
  void ICVNMapper<T, U>::beginJob()
  {}

  //......................................................................
  template <class T, class U>
  void ICVNMapper<T, U>::endJob()
  {}

  //......................................................................
  template <class T, class U>
  void ICVNMapper<T, U>::produce(art::Event& evt)
  {

    std::vector<art::Ptr<U>> hitlist;
    auto hitListHandle = evt.getHandle<std::vector<U>>(fHitsModuleLabel);
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);

    //Declaring containers for things to be stored in event
    std::unique_ptr<std::vector<cvn::PixelMap>> pmCol(new std::vector<cvn::PixelMap>);

    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    PixelMap pm = fProducer.CreateMap(detProp, hitlist);
    auto nhits = fProducer.NROI();
    pm.SetTotHits(nhits);

    if (nhits > fMinClusterHits) pmCol->push_back(pm);

    evt.put(std::move(pmCol), fClusterPMLabel);
  }

} //namespace cvn
