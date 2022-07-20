////////////////////////////////////////////////////////////////////////
// Class:       EvaluateROIEff
// Plugin Type: analyzer (art v3_05_00)
// File:        EvaluateROIEff_module.cc
//
// Generated at Sun May  3 23:16:14 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_10_00.
//
// Author: Tingjun Yang, tjyang@fnal.gov
//         Wanwei Wu, wwu@fnal.gov
//
// About: ROI efficiency and purity evaluation using "signal". Here, "signal"
//        is consecutive energy deposits based on tdc ticks.
////////////////////////////////////////////////////////////////////////

#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "TEfficiency.h"
#include "TH1D.h"

#include <cstdlib> // std::abs()
#include <string>
#include <vector>

namespace nnet {
  class EvaluateROIEff;
}

class nnet::EvaluateROIEff : public art::EDAnalyzer {
public:
  explicit EvaluateROIEff(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EvaluateROIEff(EvaluateROIEff const&) = delete;
  EvaluateROIEff(EvaluateROIEff&&) = delete;
  EvaluateROIEff& operator=(EvaluateROIEff const&) = delete;
  EvaluateROIEff& operator=(EvaluateROIEff&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  void beginJob() override;
  void endJob() override;
  bool isSignalInROI(int starttick, int endtick, int maxtick, int roistart, int roiend);
  // Declare member data here.
  art::InputTag fWireProducerLabel;
  art::InputTag
    fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector

  TH1D* h_energy[3];
  TH1D* h_energy_roi[3];

  TEfficiency* eff_energy[3];

  TH1D* h_purity[3];
  TH1D* h_purity_all;

  TH1D* h_roi[3];
  TH1D* h1_roi_max[3];
  TH1D* h1_roi_max_sim[3];

  TH1D* h1_tickdiff_max[3];

  int fCount_Roi_sig[3] = {0, 0, 0};
  int fCount_Roi_total[3] = {0, 0, 0};
};

nnet::EvaluateROIEff::EvaluateROIEff(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fWireProducerLabel(p.get<art::InputTag>("WireProducerLabel", ""))
  , fSimulationProducerLabel(p.get<art::InputTag>("SimulationProducerLabel", "largeant"))
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void
nnet::EvaluateROIEff::analyze(art::Event const& e)
{

  auto const* geo = lar::providerFrom<geo::Geometry>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);
  auto const& chStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  art::Handle<std::vector<recob::Wire>> wireListHandle;
  std::vector<art::Ptr<recob::Wire>> wires;
  if (e.getByLabel(fWireProducerLabel, wireListHandle)) {
    art::fill_ptr_vector(wires, wireListHandle);
  }

  auto simChannelHandle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

  // efficiency: according to each simulated energy deposit
  // ... Loop over simChannels
  for (auto const& channel : (*simChannelHandle)) {

    // .. get simChannel channel number
    const raw::ChannelID_t ch1 = channel.Channel();
    if (chStatus.IsBad(ch1)) continue;

    if (ch1 % 1000 == 0) mf::LogInfo("EvaluateROIEFF") << ch1;
    int view = geo->View(ch1);
    auto const& timeSlices = channel.TDCIDEMap();

    // time slice from simChannel is for individual tick
    // group neighboring time slices into a "signal"
    int tdctick_previous = -999;
    double totalE = 0.;
    double maxE = -999.;
    int maxEtick = -999;
    std::vector<int> signal_starttick;
    std::vector<int> signal_endtick;
    std::vector<double> signal_energydeposits;
    std::vector<double> signal_energy_max;
    std::vector<double> signal_max_tdctick;

    for (auto const& timeSlice : timeSlices) {
      auto const tpctime = timeSlice.first;
      auto const& energyDeposits = timeSlice.second;
      int tdctick = static_cast<int>(clockData.TPCTDC2Tick(double(tpctime)));
      if (tdctick < 0 || tdctick > int(detProp.ReadOutWindowSize()) - 1) continue;

      // for a time slice, there may exist more than one energy deposit.
      double e_deposit = 0;
      for (auto const& energyDeposit : energyDeposits) {
        e_deposit += energyDeposit.energy;
      }

      if (tdctick_previous == -999) {
        signal_starttick.push_back(tdctick);
        totalE += e_deposit;
        maxE = e_deposit;
        maxEtick = tdctick;
      }
      else if (tdctick - tdctick_previous != 1) {
        signal_starttick.push_back(tdctick);
        signal_endtick.push_back(tdctick_previous);
        signal_energydeposits.push_back(totalE);
        signal_energy_max.push_back(maxE);
        signal_max_tdctick.push_back(maxEtick);
        totalE = e_deposit;
        maxE = e_deposit;
        maxEtick = tdctick;
      }
      else if (tdctick - tdctick_previous == 1) {
        totalE += e_deposit;
        if (maxE < e_deposit) {
          maxE = e_deposit;
          maxEtick = tdctick;
        }
      }

      tdctick_previous = tdctick;

    } // loop over timeSlices timeSlice

    signal_endtick.push_back(tdctick_previous); // for last one
    signal_energydeposits.push_back(totalE);    // for last one
    signal_energy_max.push_back(maxE);          // for last one
    signal_max_tdctick.push_back(maxEtick);     // for last one

    if (signal_starttick.size() == 0 ||
        (signal_endtick.size() == 1 && signal_endtick.back() == -999))
      continue;

    for (auto& wire : wires) {
      if (wire->Channel() != ch1) continue;

      const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();

      std::vector<float> vecADC =
        wire->Signal(); // a zero-padded full length vector filled with RoI signal.

      // loop over signals:
      // a) if signal s is not in any ROI (including the case no ROI), fill h_energy 
      //    with signal_energy_max[s];
      // b) if signal s is in a ROI, check following signals that are also in this ROI,
      //    then use the maximum of signal_energy_max to fill h_energy and h_energy_roi. 
      //    After this, the loop will skip to the signal that is not in this ROI.
      for (size_t s = 0; s < signal_starttick.size(); s++) {
        // case a: signal is not in any ROI
        bool IsSignalOutside = true;
        for (const auto& range : signalROI.get_ranges()) {
          //const auto& waveform = range.data();
          raw::TDCtick_t roiFirstBinTick = range.begin_index();
          raw::TDCtick_t roiLastBinTick = range.end_index();

          if (isSignalInROI(signal_starttick[s],
                            signal_endtick[s],
                            signal_max_tdctick[s],
                            roiFirstBinTick,
                            roiLastBinTick)) {
            IsSignalOutside = false;
            break;
          }
        } // loop over range

        if (IsSignalOutside) {
          //cout << "This signal is not in any ROI: " << signal_starttick[s] << " -> " << signal_endtick[s] << endl;
          h_energy[view]->Fill(signal_energy_max[s]);
          h1_tickdiff_max[view]->Fill(-99);
          continue;
        }

        // signal is in one ROI
        for (const auto& range : signalROI.get_ranges()) {
          //const auto& waveform = range.data();
          raw::TDCtick_t roiFirstBinTick = range.begin_index();
          raw::TDCtick_t roiLastBinTick = range.end_index();

          if (isSignalInROI(signal_starttick[s],
                            signal_endtick[s],
                            signal_max_tdctick[s],
                            roiFirstBinTick,
                            roiLastBinTick)) {
            // maximum pulse height and postion for this ROI
            double maxadc_sig = 0;
            int maxadc_tick = -99;
            for (int k = roiFirstBinTick; k < roiLastBinTick; k++) {
              if (vecADC[k] > maxadc_sig) {
                maxadc_sig = vecADC[k];
                maxadc_tick = k;
              }
            }

            double maxE_roi = -999.;
            double maxE_roi_tick = -999.;

            if (maxE_roi < signal_energy_max[s]) {
              maxE_roi = signal_energy_max[s];
              maxE_roi_tick = signal_max_tdctick[s];
            }

            // check the following signals in the same ROI
            for (size_t s2 = s + 1; s2 < signal_starttick.size(); s2++) {
              if (isSignalInROI(signal_starttick[s2],
                                signal_endtick[s2],
                                signal_max_tdctick[s2],
                                roiFirstBinTick,
                                roiLastBinTick)) {
                if (maxE_roi < signal_energy_max[s2]) {
                  maxE_roi = signal_energy_max[s2];
                  maxE_roi_tick = signal_max_tdctick[s2];
                }
                if (s2 == signal_starttick.size() - 1) { s = s2; }
              }
              else {
                s = s2 - 1;
                break;
              }
            } // loop over s2

            // finish this ROI
            h_energy[view]->Fill(maxE_roi);
            h_energy_roi[view]->Fill(maxE_roi);
            h1_tickdiff_max[view]->Fill(maxE_roi_tick - maxadc_tick);
            break;
          } // isSignalInROI
        }   // loop over range

      } // loop over signals s
    }   // loop over wires wire
  }     // loop simChannels

  // purity: # signals in ROI / (#signals in ROI + #non-signals in ROI). Because we only consider the maximum signal in the ROI, this is quivalent to purity of ( number of ROIs with signal / number of ROIs)
  double roi_sig[3] = {0., 0., 0.};   // number of roi contains signal in an event
  double roi_total[3] = {0., 0., 0.}; // number of roi in an event

  for (auto& wire : wires) {
    const raw::ChannelID_t wirechannel = wire->Channel();
    if (chStatus.IsBad(wirechannel)) continue;

    int view = wire->View();

    const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();

    if (signalROI.get_ranges().empty()) continue;

    std::vector<float> vecADC =
      wire->Signal(); // a zero-padded full length vector filled with RoI signal.

    roi_total[view] += signalROI.get_ranges().size();
    fCount_Roi_total[view] += signalROI.get_ranges().size();

    for (const auto& range : signalROI.get_ranges()) {
      bool IsSigROI = false;

      raw::TDCtick_t roiFirstBinTick = range.begin_index();
      raw::TDCtick_t roiLastBinTick = range.end_index();

      double maxadc_sig = -99.;
      for (int k = roiFirstBinTick; k < roiLastBinTick; k++) {
        if (std::abs(vecADC[k]) > maxadc_sig) maxadc_sig = std::abs(vecADC[k]);
      }
      h1_roi_max[view]->Fill(maxadc_sig);

      // check simulation: ideally we could put this part outside loop over range to make the algorithm more efficient. However, as there are only one or two ROIs in a channel for most of cases, we keep this style to make the algorithm more readable. Also, signal information are kept here.
      for (auto const& channel : (*simChannelHandle)) {
        if (wirechannel != channel.Channel()) continue;

        int tdctick_previous = -999;
        double totalE = 0.;
        double maxE = -999.;
        int maxEtick = -999;
        std::vector<int> signal_starttick;
        std::vector<int> signal_endtick;
        std::vector<double> signal_energydeposits;
        std::vector<double> signal_energy_max;
        std::vector<double> signal_max_tdctick;

        auto const& timeSlices = channel.TDCIDEMap();

        for (auto const& timeSlice : timeSlices) {
          auto const tpctime = timeSlice.first;
          auto const& energyDeposits = timeSlice.second;
          int tdctick = static_cast<int>(clockData.TPCTDC2Tick(double(tpctime)));
          if (tdctick < 0 || tdctick > int(detProp.ReadOutWindowSize()) - 1) continue;

          double e_deposit = 0;
          for (auto const& energyDeposit : energyDeposits) {
            e_deposit += energyDeposit.energy;
          }

          if (tdctick_previous == -999) {
            signal_starttick.push_back(tdctick);
            totalE += e_deposit;
            maxE = e_deposit;
            maxEtick = tdctick;
          }
          else if (tdctick - tdctick_previous != 1) {
            signal_starttick.push_back(tdctick);
            signal_endtick.push_back(tdctick_previous);
            signal_energydeposits.push_back(totalE);
            signal_energy_max.push_back(maxE);
            signal_max_tdctick.push_back(maxEtick);
            totalE = e_deposit;
            maxE = e_deposit;
            maxEtick = tdctick;
          }
          else if (tdctick - tdctick_previous == 1) {
            totalE += e_deposit;
            if (maxE < e_deposit) {
              maxE = e_deposit;
              maxEtick = tdctick;
            }
          }

          tdctick_previous = tdctick;

        } // loop over timeSlices timeSlice

        signal_endtick.push_back(tdctick_previous); // for last one
        signal_energydeposits.push_back(totalE);    // for last one
        signal_energy_max.push_back(maxE);          // for last one
        signal_max_tdctick.push_back(maxEtick);     // for last one

        if (signal_starttick.size() == 0 ||
            (signal_endtick.size() == 1 && signal_endtick.back() == -999))
          continue;

        // check if signal/signals in this ROI
        for (size_t s = 0; s < signal_starttick.size(); s++) {
          if (isSignalInROI(signal_starttick[s],
                            signal_endtick[s],
                            signal_max_tdctick[s],
                            roiFirstBinTick,
                            roiLastBinTick)) {
            IsSigROI = true;
            break;
          }
        } // loop over s
      }   // loop simChannels

      h_roi[view]->Fill(0); // total roi
      if (IsSigROI) {
        roi_sig[view] += 1.;
        fCount_Roi_sig[view] += 1;
        h_roi[view]->Fill(1); // sig roi
        h1_roi_max_sim[view]->Fill(maxadc_sig);
      }
    } // loop ranges of signalROI
  }   // loop wires

  // purity of each plane for each event
  for (int i = 0; i < 3; i++) {
    if (roi_total[i]) {
      double purity = roi_sig[i] / roi_total[i];
      if (purity == 1) purity = 1. - 1.e-6;
      h_purity[i]->Fill(purity);
    }
  }

  // combined purity for each event
  if (roi_total[0] + roi_total[1] + roi_total[2]) {
    double purity =
      (roi_sig[0] + roi_sig[1] + roi_sig[2]) / (roi_total[0] + roi_total[1] + roi_total[2]);
    if (purity == 1) purity = 1. - 1.e-6;
    h_purity_all->Fill(purity);
  }
}

void
nnet::EvaluateROIEff::beginJob()
{

  art::ServiceHandle<art::TFileService const> tfs;

  for (int i = 0; i < 3; ++i) {
    h_energy[i] =
      tfs->make<TH1D>(Form("h_energy_%d", i), Form("Plane %d; Energy (MeV);", i), 100, 0, 1);
    h_energy_roi[i] =
      tfs->make<TH1D>(Form("h_energy_roi_%d", i), Form("Plane %d; Energy (MeV);", i), 100, 0, 1);

    h_purity[i] = tfs->make<TH1D>(Form("h_purity_%d", i), Form("Plane %d; Purity;", i), 20, 0, 1);
    h_roi[i] = tfs->make<TH1D>(Form("h_roi_%d", i),
                               Form("Plane %d; Purity;", i),
                               5,
                               0,
                               5); // index 0: total rois; index 1: total sig rois

    h1_roi_max[i] =
      tfs->make<TH1D>(Form("h1_roi_max_%d", i), Form("Plane %d; Max adc;", i), 50, 0, 50);
    h1_roi_max_sim[i] =
      tfs->make<TH1D>(Form("h1_roi_max_sim_%d", i), Form("Plane %d; Max adc;", i), 50, 0, 50);

    h1_tickdiff_max[i] = tfs->make<TH1D>(Form("h1_tickdiff_max_%d", i),
                                         Form("Plane %d; tick diff (maxE - max pulse);", i),
                                         100,
                                         -50,
                                         50);
  }

  h_purity_all = tfs->make<TH1D>("h_purity_all", "All Planes; Purity;", 20, 0, 1);
}

void
nnet::EvaluateROIEff::endJob()
{
  art::ServiceHandle<art::TFileService const> tfs;

  for (int i = 0; i < 3; ++i) {
    if (TEfficiency::CheckConsistency(*(h_energy_roi[i]), *(h_energy[i]))) {
      eff_energy[i] = tfs->make<TEfficiency>(*(h_energy_roi[i]), *(h_energy[i]));
      eff_energy[i]->Write(Form("eff_energy_%d", i));
    }
  }
}

bool
nnet::EvaluateROIEff::isSignalInROI(int starttick,
                                    int endtick,
                                    int maxtick,
                                    int roistart,
                                    int roiend)
{
  // For a signal in the ROI, two cases are considered:
  //   (i) signal is totally in the ROI
  //   (ii) signal is partially in the ROI
  // Equivalently, we can just check whether the maxtick is in the ROI (simple one). 
  return roistart <= maxtick && maxtick < roiend;
}

DEFINE_ART_MODULE(nnet::EvaluateROIEff)
