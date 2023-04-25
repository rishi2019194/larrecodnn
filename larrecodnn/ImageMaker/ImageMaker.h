#include "art/Framework/Principal/Event.h"
#include "hep_hpc/hdf5/File.hpp"

namespace dnn {

  class ImageMaker {
  public:
    virtual ~ImageMaker() noexcept = default;
    virtual void saveImage(art::Event const&, hep_hpc::hdf5::File&) = 0;
  };
}
