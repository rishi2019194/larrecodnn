#ifndef IWaveformDenoise_H
#define IWaveformDenoise_H

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

namespace wavdenoise_tool {
  class IWaveformDenoise {
  public:
    virtual ~IWaveformDenoise() noexcept = default;

    // Calculate multi-class probabilities for waveform
    virtual std::vector<std::vector<float>> applyDenoisingAE(
      const std::vector<std::vector<float>>&) const = 0;

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    std::vector<float> denoiseWaveform(const std::vector<float>& adcin) const
    {
      unsigned int fWaveformSize = adcin.size();
      //std::cout << " !!!! fWaveformSize: " << fWaveformSize << std::endl;

      // .. rescale input waveform for CNN
      std::vector<float> adc(fWaveformSize);
      for (size_t itck = 0; itck < fWaveformSize; ++itck) {
        adc[itck] = (adcin[itck] - fCnnMean) / fCnnScale;
      }

      float dmn =
        fWaveformSize - fWindowSize; // dist between trail edge of 1st win & last data point
      unsigned int fNumStrides =
        std::ceil(dmn / float(fStrideLength)); // # strides to scan entire waveform
      unsigned int overshoot = fNumStrides * fStrideLength + fWindowSize - fWaveformSize;
      unsigned int fLastWindowSize = fWindowSize - overshoot;

      unsigned int avgrange = fWindowSize - fStrideLength;
      unsigned int j1last, j1blast, j2last;

      if (fLastWindowOpt == 0) {       // opt0: last window shifted back so tail
        j2last = fWaveformSize;        //       edge coincides with waveformsize
        j1last = j2last - fWindowSize; //       & window is not truncated
        j1blast = (fNumStrides - 1) * fStrideLength + fWindowSize;
      }
      else { // opt1: last window is truncated
        j1last = fNumStrides * fStrideLength;
        j2last = j1last + fLastWindowSize;
        if (avgrange < fLastWindowSize) { j1blast = j1last + avgrange; }
        else {
          j1blast = j1last + fLastWindowSize;
        }
      }

      // .. create a vector of windows
      std::vector<std::vector<float>> wwv(fNumStrides + 1, std::vector<float>(fWindowSize, 0.));

      // .. fill each window with adc values
      unsigned int j1, j2, k, j1b;

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
        for (unsigned int j = j1last; j < j2last; j++) {
          wwv[fNumStrides][k] = adc[j];
          k++;
        }
      }

      // ... use 1d AE to denoise each window
      std::vector<std::vector<float>> dwv = applyDenoisingAE(wwv);

      if (fNumStrides < 1) {
        k = 0;
        for (unsigned int j = 0; j < fWindowSize; j++) {
          adc[j] = dwv[0][k];
          k++;
        }
      }
      else {
        for (unsigned int i = 0; i < fNumStrides; i++) {
          if (i == 0) {
            k = 0;
            for (unsigned int j = 0; j < fWindowSize; j++) {
              adc[j] = dwv[i][k];
              k++;
            }
          }
          else {
            j1 = i * fStrideLength;
            j1b = j1 + avgrange;
            j2 = j1 + fWindowSize;
            k = 0;
            for (unsigned int j = j1; j < j1b; j++) {
              adc[j] = 0.5 * (dwv[i][k] + adc[j]);
              k++;
            }
            for (unsigned int j = j1b; j < j2; j++) {
              adc[j] = dwv[i][k];
              k++;
            }
          }
        }
        // .. last window is a special case
        k = 0;
        for (unsigned int j = j1last; j < j1blast; j++) {
          adc[j] = 0.5 * (dwv[fNumStrides][k] + adc[j]);
          k++;
        }
        if (j1blast < j2last) {
          for (unsigned int j = j1blast; j < j2last; j++) {
            adc[j] = dwv[fNumStrides][k];
            k++;
          }
        }
      }
      // .. undo rescaling of input waveform before returning
      std::vector<float> adcout(fWaveformSize);
      for (size_t itck = 0; itck < fWaveformSize; ++itck) {
        adcout[itck] = adc[itck] * fCnnScale + fCnnMean;
      }
      return adcout;
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

    void setupWaveDenoiseParams(const fhicl::ParameterSet& pset)
    {
      // ... load the mean and scale (std) vectors
      fCnnMean = pset.get<float>("CnnMean", 0.);
      fCnnScale = pset.get<float>("CnnScale", 1.);

      fWindowSize = pset.get<unsigned int>("ScanWindowSize", 0); // 200
      fStrideLength = pset.get<unsigned int>("StrideLength", 0); // 150
      fLastWindowOpt = pset.get<unsigned int>("LastWindowOpt", 0);
    }

  private:
    float fCnnMean;
    float fCnnScale;
    unsigned int fWindowSize;   // Scan window size
    unsigned int fStrideLength; // Offset (in #time ticks) between scan windows
    unsigned int fLastWindowOpt;
  };
}

#endif
