////////////////////////////////////////////////////////////////////////
/// \file    TFNetHandler.h
/// \brief   TFNetHandler for CVN
/// \author  Leigh Whitehead
////////////////////////////////////////////////////////////////////////

#ifndef LCVN_TFNETHANDLER_H
#define LCVN_TFNETHANDLER_H

#include <memory>
#include <vector>

#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"

namespace lcvn {

  /// Wrapper for caffe::Net which handles construction and prediction
  class ITFNetHandler {
  public:

    virtual ~ITFNetHandler() noexcept = default;
    /// Return prediction arrays for PixelMap
    virtual std::vector<std::vector<float>> Predict(const PixelMap& pm) const = 0;

  };

}

#endif // CVN_TFNETHANDLER_H
