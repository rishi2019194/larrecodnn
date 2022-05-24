#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace cvn {

  typedef ICVNMapper<cvn::PixelMapHitProducer, recob::Hit> CVNMapper;
  template class ICVNMapper<cvn::PixelMapHitProducer, recob::Hit>;

DEFINE_ART_MODULE(cvn::CVNMapper)
}
