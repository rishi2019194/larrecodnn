#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"

#include "lardataobj/RecoBase/Wire.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace cvn {

  typedef ICVNMapper<cvn::PixelMapWireProducer, recob::Wire> LArCVNWireMapper;
  template class ICVNMapper<cvn::PixelMapWireProducer, recob::Wire>;

DEFINE_ART_MODULE(cvn::LArCVNWireMapper)
}
