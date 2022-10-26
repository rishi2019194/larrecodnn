#include "larrecodnn/CVN/interfaces/ICVNMapper.cxx"
#include "larrecodnn/CVN/interfaces/ICVNMapper.h"

#include "lardataobj/RecoBase/Wire.h"

#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"

namespace lcvn {

  typedef ICVNMapper<lcvn::PixelMapWireProducer, recob::Wire> LArCVNWireMapper;
  template class ICVNMapper<lcvn::PixelMapWireProducer, recob::Wire>;

  DEFINE_ART_MODULE(lcvn::LArCVNWireMapper)
}
