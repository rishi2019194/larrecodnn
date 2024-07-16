#ifndef IWaveformRecog_H
#define IWaveformRecog_H

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
    std::vector<bool> findROI(const std::vector<float>& adcin) const
    {
      std::vector<bool> bvec(fWaveformSize, false);
      if (adcin.size() != fWaveformSize) { return bvec; }

      std::vector<std::vector<float>> predv = scanWaveform(adcin);

      // .. set to true all bins in the output vector that are in windows identified as signals
      if (fNumStrides < 1) {
        if (predv[0][0] > fCnnPredCut) { std::fill_n(bvec.begin(), fWindowSize, true); }
      }
      else {
        int j1;
        for (unsigned int i = 0; i < fNumStrides; i++) {
          j1 = i * fStrideLength;
          if (predv[i][0] > fCnnPredCut) { std::fill_n(bvec.begin() + j1, fWindowSize, true); }
        }
        // .. last window is a special case
        if (predv[fNumStrides][0] > fCnnPredCut) {
          j1 = fJ1Last;
          std::fill_n(bvec.begin() + j1, fLastWindowSize, true);
        }
      }
      return bvec;
    }

    // -------------------------------------------------------------
    // Return a vector of floats of the same size as the input
    // waveform. The value in each bin represents the probability
    // whether that bin is in an ROI or not
    // -------------------------------------------------------------
    std::vector<float> predROI(const std::vector<float>& adcin) const
    {
      std::vector<float> fvec(fWaveformSize, 0.);
      if (adcin.size() != fWaveformSize) { return fvec; }

      std::vector<std::vector<float>> predv = scanWaveform(adcin);

      // .. set value in each bin of output vector to the prediction for the window it is in
      if (fNumStrides < 1) { std::fill_n(fvec.begin(), fWindowSize, predv[0][0]); }
      else {
        int j1;
        for (unsigned int i = 0; i < fNumStrides; i++) {
          j1 = i * fStrideLength;
          std::fill_n(fvec.begin() + j1, fWindowSize, predv[i][0]);
        }
        // .. last window is a special case
        j1 = fJ1Last;
        std::fill_n(fvec.begin() + j1, fLastWindowSize, predv[fNumStrides][0]);
      }
      return fvec;
    }

    // -------------------------------------------------------------
    // Like findROI but accepts any length input waveform
    // -------------------------------------------------------------
    std::vector<bool> findROI2(const std::vector<float>& adcin) const
    {
      unsigned int waveformsize = adcin.size();

      // .. rescale input waveform for CNN
      std::vector<float> adc(waveformsize);
      for (size_t itck = 0; itck < waveformsize; ++itck) {
        adc[itck] = (adcin[itck] - fCnnMean) / fCnnScale;
      }

      float dmn =
        waveformsize - fWindowSize; // dist between trail edge of 1st win & last data point
      unsigned int numstrides =
        std::ceil(dmn / float(fStrideLength)); // # strides to scan entire waveform
      unsigned int overshoot = numstrides * fStrideLength + fWindowSize - waveformsize;
      unsigned int lastwindowsize = fWindowSize - overshoot;

      unsigned int j1last, j2last;

      if (fLastWindowOpt == 0) {       // opt0: last window shifted back so tail
        j2last = waveformsize;         //       edge coincides with waveformsize
        j1last = j2last - fWindowSize; //       & window is not truncated
        lastwindowsize = fWindowSize;
      }
      else { // opt1: last window is truncated
        j1last = numstrides * fStrideLength;
        j2last = j1last + lastwindowsize;
      }

      // .. create a vector of windows
      std::vector<std::vector<float>> wwv(numstrides + 1, std::vector<float>(fWindowSize, 0.));

      // .. fill each window with adc values
      unsigned int j1, j2, k;
      if (numstrides < 1) {
        k = 0;
        for (unsigned int j = 0; j < fWindowSize; j++) {
          wwv[0][k] = adc[j];
          k++;
        }
      }
      else {
        for (unsigned int i = 0; i < numstrides; i++) {
          j1 = i * fStrideLength;
          j2 = j1 + fWindowSize;
          k = 0;
          for (unsigned int j = j1; j < j2; j++) {
            wwv[i][k] = adc[j];
            k++;
          }
        }
        // .. last window is a special case
        k = 0;
        for (unsigned int j = j1last; j < j2last; j++) {
          wwv[numstrides][k] = adc[j];
          k++;
        }
      }
      // ... use waveform recognition CNN to perform inference on each window
      std::vector<std::vector<float>> predv = predictWaveformType(wwv);

      // .. set to true all bins in the output vector that are in windows identified as signals
      std::vector<bool> bvec(waveformsize, false);
      if (numstrides < 1) {
        if (predv[0][0] > fCnnPredCut) { std::fill_n(bvec.begin(), fWindowSize, true); }
      }
      else {
        int jj;
        for (unsigned int i = 0; i < numstrides; i++) {
          jj = i * fStrideLength;
          if (predv[i][0] > fCnnPredCut) { std::fill_n(bvec.begin() + jj, fWindowSize, true); }
        }
        // .. last window is a special case
        if (predv[numstrides][0] > fCnnPredCut) {
          jj = j1last;
          std::fill_n(bvec.begin() + jj, lastwindowsize, true);
        }
      }
      return bvec;
    }

    void configureWaveformSize(unsigned int waveformsize)
    {
      fWaveformSize = waveformsize;
      if (fWaveformSize > 0 && fWindowSize > 0) {
        float dmn =
          fWaveformSize - fWindowSize; // dist between trail edge of 1st win & last data point
        fNumStrides = std::ceil(dmn / float(fStrideLength)); // # strides to scan entire waveform
        unsigned int overshoot = fNumStrides * fStrideLength + fWindowSize - fWaveformSize;
        fLastWindowSize = fWindowSize - overshoot;
        //unsigned int numwindows = fNumStrides + 1;
        if (fLastWindowOpt == 0) {         // opt0: last window shifted back so tail
          fJ2Last = fWaveformSize;         //       edge coincides with waveformsize
          fJ1Last = fJ2Last - fWindowSize; //       & window is not truncated
          fLastWindowSize = fWindowSize;
        }
        else { // opt1: last window is truncated
          fJ1Last = fNumStrides * fStrideLength;
          fJ2Last = fJ1Last + fLastWindowSize;
        }
        //std::cout << " !!!!! WaveformRecog: WaveformSize = " << fWaveformSize
        //          << ", WindowSize = " << fWindowSize
        //          << ", StrideLength = " << fStrideLength
        //          << ", dmn/StrideLength = " << dmn / fStrideLength << std::endl;
        //std::cout << "       dmn = " << dmn << ", NumStrides = " << fNumStrides
        //          << ", overshoot = " << overshoot << ", LastWindowSize = " << fLastWindowSize
        //          << ", numwindows = " << numwindows << std::endl;
      }
    }

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

    void setupWaveRecRoiParams(const fhicl::ParameterSet& pset)
    {
      fCnnPredCut = pset.get<float>("CnnPredCut", 0.5);
      unsigned int waveform_size = pset.get<unsigned int>("WaveformSize", 0); // 6000

      // ... load the mean and scale (std) vectors
      fCnnMean = pset.get<float>("CnnMean", 0.);
      fCnnScale = pset.get<float>("CnnScale", 1.);

      fWindowSize = pset.get<unsigned int>("ScanWindowSize", 0); // 200
      fStrideLength = pset.get<unsigned int>("StrideLength", 0); // 150
      fLastWindowOpt = pset.get<unsigned int>("LastWindowOpt", 0);

      configureWaveformSize(waveform_size);
    }

  private:
    float fCnnMean;
    float fCnnScale;
    float fCnnPredCut;
    unsigned int fWaveformSize; // Full waveform size
    unsigned int fWindowSize;   // Scan window size
    unsigned int fStrideLength; // Offset (in #time ticks) between scan windows
    unsigned int fLastWindowOpt;
    unsigned int fNumStrides;
    unsigned int fLastWindowSize;
    unsigned int fJ1Last;
    unsigned int fJ2Last;

    std::vector<std::vector<float>> scanWaveform(const std::vector<float>& adcin) const
    {
      // .. rescale input waveform for CNN
      std::vector<float> adc(fWaveformSize);
      for (size_t itck = 0; itck < fWaveformSize; ++itck) {
        adc[itck] = (adcin[itck] - fCnnMean) / fCnnScale;
      }

      // .. create a vector of windows
      std::vector<std::vector<float>> wwv(fNumStrides + 1, std::vector<float>(fWindowSize, 0.));

      // .. fill each window with adc values
      unsigned int j1, j2, k;
      if (fNumStrides < 1) {
        k = 0;
        for (unsigned int j = 0; j < fWindowSize; j++) {
          wwv[0][k] = adc[j];
          k++;
        }
      }
      else {
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
        k = 0;
        for (unsigned int j = fJ1Last; j < fJ2Last; j++) {
          wwv[fNumStrides][k] = adc[j];
          k++;
        }
      }

      // ... use waveform recognition CNN to perform inference on each window
      return predictWaveformType(wwv);
    }
  };
}

#endif
