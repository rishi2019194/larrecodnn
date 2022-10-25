#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"

#include "lardataobj/RecoBase/Hit.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace cvn {

  typedef ICVNMapper<cvn::PixelMapHitProducer, recob::Hit> LArCVNHitMapper;
  template class ICVNMapper<cvn::PixelMapHitProducer, recob::Hit>;

  DEFINE_ART_MODULE(cvn::LArCVNHitMapper)
}
