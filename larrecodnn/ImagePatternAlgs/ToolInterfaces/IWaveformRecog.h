#ifndef IWaveformRecog_H
#define IWaveformRecog_H

#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace wavrec_tool {
  class IWaveformRecog {
  public:
    virtual ~IWaveformRecog() noexcept = default;

    // Calculate multi-class probabilities for waveform
    virtual std::vector<std::vector<float>> predictWaveformType(
      const std::vector<std::vector<float>>&) const = 0;

    // ---------------------------------------------------------------------
    // Return a vector of booleans of the same size as the input  waveform.
    // The value of each element of the vector represents whether the
    // corresponding time bin of the waveform is in an ROI or not.
    // ---------------------------------------------------------------------
    std::vector<bool>
    findROI(const std::vector<float>& adcin) const
    {
      std::vector<bool> bvec(fWaveformSize, false);
      if (adcin.size() != fWaveformSize) { return bvec; }

      std::vector<std::vector<float>> predv = scanWaveform(adcin);

      // .. set to true all bins in the output vector that are in windows identified as signals
      int j1;
      for (unsigned int i = 0; i < fNumStrides; i++) {
        j1 = i * fStrideLength;
        if (predv[i][0] > fCnnPredCut) { std::fill_n(bvec.begin() + j1, fWindowSize, true); }
      }
      // .. last window is a special case
      if (predv[fNumStrides][0] > fCnnPredCut) {
        j1 = fNumStrides * fStrideLength;
        std::fill_n(bvec.begin() + j1, fLastWindowSize, true);
      }
      return bvec;
    }

    // -------------------------------------------------------------
    // Return a vector of floats of the same size as the input
    // waveform. The value in each bin represents the probability
    // whether that bin is in an ROI or not
    // -------------------------------------------------------------
    std::vector<float>
    predROI(const std::vector<float>& adcin) const
    {
      std::vector<float> fvec(fWaveformSize, 0.);
      if (adcin.size() != fWaveformSize) { return fvec; }

      std::vector<std::vector<float>> predv = scanWaveform(adcin);

      // .. set value in each bin of output vector to the prediction for the window it is in
      int j1;
      for (unsigned int i = 0; i < fNumStrides; i++) {
        j1 = i * fStrideLength;
        std::fill_n(fvec.begin() + j1, fWindowSize, predv[i][0]);
      }
      // .. last window is a special case
      j1 = fNumStrides * fStrideLength;
      std::fill_n(fvec.begin() + j1, fLastWindowSize, predv[fNumStrides][0]);
      return fvec;
    }

  protected:
    std::string
    findFile(const char* fileName) const
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

    void
    setupWaveRecRoiParams(const fhicl::ParameterSet& pset)
    {
      fCnnPredCut = pset.get<float>("CnnPredCut", 0.5);
      fWaveformSize = pset.get<unsigned int>("WaveformSize", 0); // 6000
      std::string fMeanFilename = pset.get<std::string>("MeanFilename", "");
      std::string fScaleFilename = pset.get<std::string>("ScaleFilename", "");

      // ... load the mean and scale (std) vectors
      if (!fMeanFilename.empty() && !fScaleFilename.empty()) {
        float val;
        std::ifstream meanfile(findFile(fMeanFilename.c_str()).c_str());
        if (meanfile.is_open()) {
          while (meanfile >> val)
            meanvec.push_back(val);
          meanfile.close();
          if (meanvec.size() != fWaveformSize) {
            throw cet::exception("WaveformRecogTf_tool")
              << "vector of mean values does not match waveform size, exiting" << std::endl;
          }
        }
        else {
          throw cet::exception("WaveformRecogTf_tool")
            << "failed opening StdScaler mean file, exiting" << std::endl;
        }
        std::ifstream scalefile(findFile(fScaleFilename.c_str()).c_str());
        if (scalefile.is_open()) {
          while (scalefile >> val)
            scalevec.push_back(val);
          scalefile.close();
          if (scalevec.size() != fWaveformSize) {
            throw cet::exception("WaveformRecogTf_tool")
              << "vector of scale values does not match waveform size, exiting" << std::endl;
          }
        }
        else {
          throw cet::exception("WaveformRecogTf_tool")
            << "failed opening StdScaler scale file, exiting" << std::endl;
        }
      }
      else {
        fCnnMean = pset.get<float>("CnnMean", 0.);
        fCnnScale = pset.get<float>("CnnScale", 1.);
        meanvec.resize(fWaveformSize);
        std::fill(meanvec.begin(), meanvec.end(), fCnnMean);
        scalevec.resize(fWaveformSize);
        std::fill(scalevec.begin(), scalevec.end(), fCnnScale);
      }

      fWindowSize = pset.get<unsigned int>("ScanWindowSize", 0); // 200
      fStrideLength = pset.get<unsigned int>("StrideLength", 0); // 150

      if (fWaveformSize > 0 && fWindowSize > 0) {
        float dmn =
          fWaveformSize - fWindowSize; // dist between trail edge of 1st win & last data point
        fNumStrides = std::ceil(dmn / float(fStrideLength)); // # strides to scan entire waveform
        unsigned int overshoot = fNumStrides * fStrideLength + fWindowSize - fWaveformSize;
        fLastWindowSize = fWindowSize - overshoot;
        unsigned int numwindows = fNumStrides + 1;
        std::cout << " !!!!! WaveformRoiFinder: WindowSize = " << fWindowSize
                  << ", StrideLength = " << fStrideLength
                  << ", dmn/StrideLength = " << dmn / fStrideLength << std::endl;
        std::cout << "	 dmn = " << dmn << ", NumStrides = " << fNumStrides
                  << ", overshoot = " << overshoot << ", LastWindowSize = " << fLastWindowSize
                  << ", numwindows = " << numwindows << std::endl;
      }
    }

  private:
    std::vector<float> scalevec;
    std::vector<float> meanvec;
    float fCnnMean;
    float fCnnScale;
    float fCnnPredCut;
    unsigned int fWaveformSize; // Full waveform size
    unsigned int fWindowSize;   // Scan window size
    unsigned int fStrideLength; // Offset (in #time ticks) between scan windows
    unsigned int fNumStrides;
    unsigned int fLastWindowSize;

    std::vector<std::vector<float>>
    scanWaveform(const std::vector<float>& adcin) const
    {
      // .. rescale input waveform for CNN
      std::vector<float> adc(fWaveformSize);
      for (size_t itck = 0; itck < fWaveformSize; ++itck) {
        adc[itck] = (adcin[itck] - meanvec[itck]) / scalevec[itck];
      }

      // .. create a vector of windows
      std::vector<std::vector<float>> wwv(fNumStrides + 1, std::vector<float>(fWindowSize, 0.));

      // .. fill each window with adc values
      unsigned int j1, j2, k;
      for (unsigned int i = 0; i < fNumStrides; i++) {
        j1 = i * fStrideLength;
        j2 = j1 + fWindowSize;
        k = 0;
        for (unsigned int j = j1; j < j2; j++) {
          wwv[i][k] = adc[j];
          k++;
        }
      }
      // .. last window is a special case
      j1 = fNumStrides * fStrideLength;
      j2 = j1 + fLastWindowSize;
      k = 0;
      for (unsigned int j = j1; j < j2; j++) {
        wwv[fNumStrides][k] = adc[j];
        k++;
      }

      // ... use waveform recognition CNN to perform inference on each window
      return predictWaveformType(wwv);
    }
  };
}

#endif
