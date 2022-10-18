#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"

#include "lardataobj/Simulation/SimChannel.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace cvn {

  typedef ICVNMapper<cvn::PixelMapSimProducer, sim::SimChannel> LArCVNSimMapper;
  template class ICVNMapper<cvn::PixelMapSimProducer, sim::SimChannel>;

DEFINE_ART_MODULE(cvn::LArCVNSimMapper)
}
