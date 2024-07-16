#ifndef IWireframeRecog_H
#define IWireframeRecog_H

#include "canvas/Utilities/Exception.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace wframerec_tool {
  class IWireframeRecog {
  public:
    virtual ~IWireframeRecog() noexcept = default;

    // Calculate multi-class probabilities for raw wireframes
    virtual std::vector<std::vector<float>> predictWireframeType(
      const std::vector<std::vector<std::vector<short>>>&) const = 0;

    // ---------------------------------------------------------------------
    // Return a vector of booleans of the same size as the full input frame.
    // The value of each element of the vector represents whether the
    // corresponding time bin of the full input frame is in an ROI or not.
    // ---------------------------------------------------------------------
    std::vector<bool> findROI(const std::vector<std::vector<short>>& adcin) const
    {
      std::vector<bool> bvec(fFullReadoutDepth, false);
      if (adcin.front().size() != fFullReadoutDepth) { return bvec; }

      std::vector<std::vector<float>> predv = scanFullReadout(adcin);

      // .. set to true all bins in the output vector that are in windows identified as signals
      int j1;
      for (unsigned int i = 0; i < fNumStrides; i++) {
        j1 = i * fStrideLength;
        if (predv[i][1] > fCnnPredCut) { std::fill_n(bvec.begin() + j1, fFrameWidth, true); }
      }
      // .. last window is a special case
      if (predv[fNumStrides][1] > fCnnPredCut) {
        j1 = fNumStrides * fStrideLength;
        std::fill_n(bvec.begin() + j1, fLastFrameWidth, true);
      }
      return bvec;
    }

    // -------------------------------------------------------------
    // Return a vector of floats of the same size as the full input
    // wireframe. The value in each bin represents the probability
    // whether that bin is in an ROI or not
    // -------------------------------------------------------------
    std::vector<float> predROI(const std::vector<std::vector<short>>& adcin) const
    {
      std::vector<float> fvec(fFullReadoutDepth, 0.);
      if (adcin.front().size() != fFullReadoutDepth) { return fvec; }

      std::vector<std::vector<float>> predv = scanFullReadout(adcin);

      // .. set value in each bin of output vector to the prediction for the window it is in
      int j1;
      for (unsigned int i = 0; i < fNumStrides; i++) {
        j1 = i * fStrideLength;
        std::fill_n(fvec.begin() + j1, fFrameWidth, predv[i][1]);
      }
      // .. last window is a special case
      j1 = fNumStrides * fStrideLength;
      std::fill_n(fvec.begin() + j1, fLastFrameWidth, predv[fNumStrides][1]);
      return fvec;
    }

    int countPassingFrames(const std::vector<std::vector<short>>& adcin) const
    {
      if (adcin.front().size() != fFullReadoutDepth) { return -1; }

      int npass = 0;
      std::vector<std::vector<float>> predv = scanFullReadout(adcin);
      for (unsigned int i = 0; i <= fNumStrides; i++) {
        if (predv[i][1] > fCnnPredCut) { npass++; }
      }
      return npass;
    }

    int countFailingFrames(const std::vector<std::vector<short>>& adcin) const
    {
      if (adcin.front().size() != fFullReadoutDepth) { return -1; }

      int nfail = 0;
      std::vector<std::vector<float>> predv = scanFullReadout(adcin);
      for (unsigned int i = 0; i <= fNumStrides; i++) {
        if (predv[i][1] <= fCnnPredCut) { nfail++; }
      }
      return nfail;
    }

    float GetCnnPredCut() const { return fCnnPredCut; }
    unsigned int GetScanFrameWidth() const { return fFrameWidth; }
    unsigned int GetStrideLength() const { return fStrideLength; }
    unsigned int GetNumStrides() const { return fNumStrides; }
    unsigned int GetLastFrameWidth() const { return fLastFrameWidth; }

  protected:
    std::string findFile(const char* fileName) const
    {
      std::string fname_out;
      cet::search_path sp("FW_SEARCH_PATH");
      if (!sp.find_file(fileName, fname_out)) {
        struct stat buffer;
        if (stat(fileName, &buffer) == 0) { fname_out = fileName; }
        else {
          throw art::Exception(art::errors::NotFound)
            << "Could not find the model file " << fileName;
        }
      }
      return fname_out;
    }

    void setupWframeRecRoiParams(const fhicl::ParameterSet& pset)
    {
      fCnnPredCut = pset.get<float>("CnnPredCut", 0.5);
      fFullReadoutDepth = pset.get<unsigned int>("FullReadoutDepth", 0); // 6000

      fFrameWidth = pset.get<unsigned int>("ScanFrameWidth", 0); // 200
      fStrideLength = pset.get<unsigned int>("StrideLength", 0); // 150
      fLastWindowOpt = pset.get<unsigned int>("LastWindowOpt", 0);

      if (fFullReadoutDepth > 0 && fFrameWidth > 0) {
        float dmn =
          fFullReadoutDepth - fFrameWidth; // dist between trail edge of 1st win & last data point
        fNumStrides = std::ceil(dmn / float(fStrideLength)); // # strides to scan entire waveform
        unsigned int overshoot = fNumStrides * fStrideLength + fFrameWidth - fFullReadoutDepth;
        fLastFrameWidth = fFrameWidth - overshoot;
        unsigned int numwindows = fNumStrides + 1;

        if (fLastWindowOpt == 0) {         // opt0: last window shifted back so tail
          fJ2Last = fFullReadoutDepth;     //       edge coincides with waveformsize
          fJ1Last = fJ2Last - fFrameWidth; //       & window is not truncated
        }
        else { // opt1: last window is truncated
          fJ1Last = fNumStrides * fStrideLength;
          fJ2Last = fJ1Last + fLastFrameWidth;
        }
        std::cout << " !!!!! WireframeRecog: FullReadoutDepth = " << fFullReadoutDepth
                  << ", FrameWidth = " << fFrameWidth << ", StrideLength = " << fStrideLength
                  << ", dmn/StrideLength = " << dmn / fStrideLength << std::endl;
        std::cout << "	 dmn = " << dmn << ", NumStrides = " << fNumStrides
                  << ", overshoot = " << overshoot << ", LastFrameWidth = " << fLastFrameWidth
                  << ", numwindows = " << numwindows << std::endl;
        std::cout << ", LastWindowOpt = " << fLastWindowOpt << ", Last frame indices = [" << fJ1Last
                  << "," << fJ2Last << "]" << std::endl;
      }
    }

  private:
    float fCnnPredCut;
    unsigned int fFullReadoutDepth; // Full waveform size
    unsigned int fFrameWidth;       // Scan window size
    unsigned int fStrideLength;     // Offset (in #time ticks) between scan windows
    unsigned int fNumStrides;
    unsigned int fLastFrameWidth;
    unsigned int fLastWindowOpt;
    unsigned int fJ1Last, fJ2Last;

    std::vector<std::vector<float>> scanFullReadout(
      const std::vector<std::vector<short>>& adcin) const
    {
      // .. create a vector of frames
      unsigned int rows = adcin.size();
      std::vector<std::vector<std::vector<short>>> wwv(
        fNumStrides + 1,
        std::vector<std::vector<short>>(rows, std::vector<short>(fFrameWidth, 0.)));

      // .. fill each window with adc values
      unsigned int j1, j2, k;
      for (unsigned int i = 0; i < fNumStrides; i++) {
        j1 = i * fStrideLength;
        j2 = j1 + fFrameWidth;
        for (unsigned int r = 0; r < rows; r++) {
          k = 0;
          for (unsigned int j = j1; j < j2; j++) {
            wwv[i][r][k] = adcin[r][j];
            k++;
          }
        }
      }
      // .. last window is a special case
      for (unsigned int r = 0; r < rows; r++) {
        k = 0;
        for (unsigned int j = fJ1Last; j < fJ2Last; j++) {
          wwv[fNumStrides][r][k] = adcin[r][j];
          k++;
        }
      }

      // ... use waveform recognition CNN to perform inference on each window
      return predictWireframeType(wwv);
    }
  };
}

#endif
