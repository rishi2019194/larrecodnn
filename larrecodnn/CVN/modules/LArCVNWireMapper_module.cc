#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"

#include "lardataobj/RecoBase/Wire.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace lcvn {

  typedef ICVNMapper<lcvn::PixelMapWireProducer, recob::Wire> LArCVNWireMapper;
  template class ICVNMapper<lcvn::PixelMapWireProducer, recob::Wire>;

  DEFINE_ART_MODULE(lcvn::LArCVNWireMapper)
}
