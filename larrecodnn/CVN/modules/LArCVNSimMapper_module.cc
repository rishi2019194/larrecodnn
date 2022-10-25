#include "larrecodnn/CVN/module_helpers/ICVNMapper.cxx"
#include "larrecodnn/CVN/module_helpers/ICVNMapper.h"

#include "lardataobj/Simulation/SimChannel.h"

#include "larrecodnn/CVN/module_helpers/PixelMapProducer.h"

namespace lcvn {

  typedef ICVNMapper<lcvn::PixelMapSimProducer, sim::SimChannel> LArCVNSimMapper;
  template class ICVNMapper<lcvn::PixelMapSimProducer, sim::SimChannel>;

  DEFINE_ART_MODULE(lcvn::LArCVNSimMapper)
}
