#include "larrecodnn/CVN/interfaces/ICVNMapper.cxx"
#include "larrecodnn/CVN/interfaces/ICVNMapper.h"

#include "lardataobj/RecoBase/Hit.h"

#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"

namespace lcvn {

  typedef ICVNMapper<lcvn::PixelMapHitProducer, recob::Hit> LArCVNHitMapper;
  template class ICVNMapper<lcvn::PixelMapHitProducer, recob::Hit>;

  DEFINE_ART_MODULE(lcvn::LArCVNHitMapper)
}
