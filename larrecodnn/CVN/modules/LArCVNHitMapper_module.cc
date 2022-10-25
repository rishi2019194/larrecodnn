#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"

#include "lardataobj/RecoBase/Hit.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace lcvn {

  typedef ICVNMapper<lcvn::PixelMapHitProducer, recob::Hit> LArCVNHitMapper;
  template class ICVNMapper<lcvn::PixelMapHitProducer, recob::Hit>;

  DEFINE_ART_MODULE(lcvn::LArCVNHitMapper)
}
