#include "larrecodnn/CVN/interfaces/ICVNMapper.cxx"
#include "larrecodnn/CVN/interfaces/ICVNMapper.h"

#include "lardataobj/Simulation/SimChannel.h"

#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"

namespace lcvn {

  typedef ICVNMapper<lcvn::PixelMapSimProducer, sim::SimChannel> LArCVNSimMapper;
  template class ICVNMapper<lcvn::PixelMapSimProducer, sim::SimChannel>;

  DEFINE_ART_MODULE(lcvn::LArCVNSimMapper)
}
