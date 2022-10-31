////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdTrainingData
// Author:      P.Plonski, R.Sulej (Robert.Sulej@cern.ch), D.Stefan, May 2016
//
// Training data for PointIdAlg
//
//      We use this to dump deconv. ADC for preparation of various classifiers.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larrecodnn/ImagePatternAlgs/Modules/c2numpy.h"
#include "larrecodnn/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"

#include "TH2F.h" // ADC and deposit maps
#include "TH2I.h" // PDG+vertex info map

// C++ Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
  template <typename Hist>
  void writeAndDelete(Hist*& hist)
  {
    if (!hist) return;
    hist->Write();
    delete hist;
    hist = nullptr;
  } // writeAndDelete()
} // local namespace

namespace nnet {

  class PointIdTrainingData : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<nnet::TrainingDataAlg::Config> TrainingDataAlg{Name("TrainingDataAlg")};

      fhicl::Atom<std::string> OutTextFilePath{Name("OutTextFilePath"),
                                               Comment("Text files with all needed data dumped.")};

      fhicl::Atom<std::string> OutNumpyFileName{Name("OutNumpyFileName"),
                                                Comment("Numpy files with patches.")};

      fhicl::Atom<bool> DumpToRoot{
        Name("DumpToRoot"),
        Comment("Dump to ROOT histogram file (replaces the text files)")};

      fhicl::Atom<bool> DumpToNumpy{Name("DumpToNumpy"),
                                    Comment("Dump to Numpy file (replaces the text files)")};

      fhicl::Sequence<int> SelectedTPC{
        Name("SelectedTPC"),
        Comment("use selected views only, or all views if empty list")};

      fhicl::Sequence<int> SelectedView{
        Name("SelectedView"),
        Comment("use selected tpc's only, or all tpc's if empty list")};

      fhicl::Atom<bool> Crop{Name("Crop"),
                             Comment("Crop the projection to the event region plus margin")};

      fhicl::Atom<int> Patch_size_w{Name("Patch_size_w"), Comment("Patch size in wire dimension")};

      fhicl::Atom<int> Patch_size_d{Name("Patch_size_d"), Comment("Patch size in drift dimension")};

      fhicl::Atom<double> Em{Name("Em"), Comment("Fraction of Em patches to keep")};

      fhicl::Atom<double> Trk{Name("Trk"), Comment("Fraction of Trk patches to keep")};

      fhicl::Atom<double> Michel{Name("Michel"), Comment("Fraction of Michel patches to keep")};

      fhicl::Atom<double> None{Name("None"), Comment("Fraction of None patches to keep")};

      fhicl::Atom<double> StopTrk{Name("StopTrk"),
                                  Comment("Fraction of stopping Trk patches to keep")};

      fhicl::Atom<double> CleanTrk{Name("CleanTrk"),
                                   Comment("Fraction of clean track patches to keep")};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit PointIdTrainingData(Parameters const& config);

    void beginJob() override;
    void endJob() override;

  private:
    void analyze(const art::Event& event) override;

    nnet::TrainingDataAlg fTrainingDataAlg;

    std::string fOutTextFilePath;
    std::string fOutNumpyFileName;
    bool fDumpToRoot;
    bool fDumpToNumpy;

    std::vector<int> fSelectedTPC;
    std::vector<int> fSelectedPlane;

    int fEvent;  /// number of the event being processed
    int fRun;    /// number of the run being processed
    int fSubRun; /// number of the sub-run being processed

    bool fCrop; /// crop data to event (set to false when dumping noise!)

    int fPatch_size_w; /// patch size in wire dimension
    int fPatch_size_d; /// patch size in drift dimension

    double fEm, fTrk, fMichel, fNone;
    double fStopTrk, fCleanTrk;

    int nEm, nTrk, nMichel, nNone;
    int nEm_sel, nTrk_sel, nMichel_sel, nNone_sel;
    int nStopTrk_sel, nCleanTrk_sel;

    geo::GeometryCore const* fGeometry;

    c2numpy_writer npywriter;

    CLHEP::HepRandomEngine& fEngine; ///< art-managed random-number engine

    int WeightedFit(int n,
                    std::vector<double> const& x,
                    std::vector<double> const& y,
                    std::vector<double> const& w,
                    double* params) const;
  };

  //-----------------------------------------------------------------------
  PointIdTrainingData::PointIdTrainingData(PointIdTrainingData::Parameters const& config)
    : art::EDAnalyzer(config)
    , fTrainingDataAlg(config().TrainingDataAlg())
    , fOutTextFilePath(config().OutTextFilePath())
    , fOutNumpyFileName(config().OutNumpyFileName())
    , fDumpToRoot(config().DumpToRoot())
    , fDumpToNumpy(config().DumpToNumpy())
    , fSelectedTPC(config().SelectedTPC())
    , fSelectedPlane(config().SelectedView())
    , fCrop(config().Crop())
    , fPatch_size_w(config().Patch_size_w())
    , fPatch_size_d(config().Patch_size_d())
    , fEm(config().Em())
    , fTrk(config().Trk())
    , fMichel(config().Michel())
    , fNone(config().None())
    , fStopTrk(config().StopTrk())
    , fCleanTrk(config().CleanTrk())
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this))
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry const>());

    if (fSelectedTPC.empty()) {
      for (size_t tpc = 0; tpc < fGeometry->NTPC(); ++tpc)
        fSelectedTPC.push_back(tpc);
    }

    if (fSelectedPlane.empty()) {
      for (size_t p = 0; p < fGeometry->MaxPlanes(); ++p)
        fSelectedPlane.push_back(p);
    }
  }

  //-----------------------------------------------------------------------
  void PointIdTrainingData::beginJob()
  {

    c2numpy_init(&npywriter, fOutNumpyFileName, 50000);
    c2numpy_addcolumn(&npywriter, "run", C2NUMPY_UINT32);
    c2numpy_addcolumn(&npywriter, "subrun", C2NUMPY_UINT32);
    c2numpy_addcolumn(&npywriter, "evt", C2NUMPY_UINT32);
    c2numpy_addcolumn(&npywriter, "tpc", C2NUMPY_UINT8);
    c2numpy_addcolumn(&npywriter, "plane", C2NUMPY_UINT8);
    c2numpy_addcolumn(&npywriter, "wire", C2NUMPY_UINT16);
    c2numpy_addcolumn(&npywriter, "tck", C2NUMPY_UINT16);
    c2numpy_addcolumn(&npywriter, "y0", C2NUMPY_UINT8);
    c2numpy_addcolumn(&npywriter, "y1", C2NUMPY_UINT8);
    c2numpy_addcolumn(&npywriter, "y2", C2NUMPY_UINT8);
    c2numpy_addcolumn(&npywriter, "y3", C2NUMPY_UINT8);
    for (int i = 0; i < fPatch_size_w * fPatch_size_d; ++i) {
      c2numpy_addcolumn(&npywriter, Form("x%d", i), C2NUMPY_FLOAT32);
    }
    nEm = 0;
    nTrk = 0;
    nMichel = 0;
    nNone = 0;
    nEm_sel = 0;
    nTrk_sel = 0;
    nMichel_sel = 0;
    nNone_sel = 0;
    nStopTrk_sel = 0;
    nCleanTrk_sel = 0;
  }

  //-----------------------------------------------------------------------
  void PointIdTrainingData::endJob()
  {
    std::cout << "nEm = " << nEm << std::endl;
    std::cout << "nTrk = " << nTrk << std::endl;
    std::cout << "nMichel = " << nMichel << std::endl;
    std::cout << "nNone = " << nNone << std::endl;
    std::cout << std::endl;
    std::cout << "nEm_sel = " << nEm_sel << std::endl;
    std::cout << "nTrk_sel = " << nTrk_sel << std::endl;
    std::cout << "nMichel_sel = " << nMichel_sel << std::endl;
    std::cout << "nNone_sel = " << nNone_sel << std::endl;
    std::cout << "nStopTrk_sel = " << nStopTrk_sel << std::endl;
    std::cout << "nCleanTrk_sel = " << nCleanTrk_sel << std::endl;
    c2numpy_close(&npywriter);
  }

  //-----------------------------------------------------------------------
  void PointIdTrainingData::analyze(const art::Event& event)
  {
    fEvent = event.id().event();
    fRun = event.run();
    fSubRun = event.subRun();

    bool saveSim = fTrainingDataAlg.saveSimInfo() && !event.isRealData();

    std::ostringstream os;
    os << "event_" << fEvent << "_run_" << fRun << "_subrun_" << fSubRun;

    std::cout << "analyze " << os.str() << std::endl;

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

    CLHEP::RandFlat flat(fEngine);

    for (size_t i = 0; i < fSelectedTPC.size(); ++i)
      for (size_t v = 0; v < fSelectedPlane.size(); ++v) {
        fTrainingDataAlg.setEventData(
          event, clockData, detProp, fSelectedPlane[v], fSelectedTPC[i], 0);

        unsigned int w0, w1, d0, d1;
        if (fCrop && saveSim) {
          if (fTrainingDataAlg.findCrop(0.004F, w0, w1, d0, d1)) {
            std::cout << "   crop: " << w0 << " " << w1 << " " << d0 << " " << d1 << std::endl;
          }
          else {
            std::cout << "   skip empty tpc:" << fSelectedTPC[i] << " / view:" << fSelectedPlane[v]
                      << std::endl;
            continue;
          }
        }
        else {
          w0 = 0;
          w1 = fTrainingDataAlg.NWires();
          d0 = 0;
          d1 = fTrainingDataAlg.NScaledDrifts();
        }

        if (fDumpToRoot) {
          std::ostringstream ss1;
          ss1 << "raw_" << os.str() << "_tpc_" << fSelectedTPC[i] << "_view_"
              << fSelectedPlane[v]; // TH2's name

          art::ServiceHandle<art::TFileService const> tfs;
          TH2F* rawHist =
            tfs->make<TH2F>((ss1.str() + "_raw").c_str(), "ADC", w1 - w0, w0, w1, d1 - d0, d0, d1);
          TH2F* depHist = 0;
          TH2I* pdgHist = 0;
          if (saveSim) {
            depHist = tfs->make<TH2F>(
              (ss1.str() + "_deposit").c_str(), "Deposit", w1 - w0, w0, w1, d1 - d0, d0, d1);
            pdgHist = tfs->make<TH2I>(
              (ss1.str() + "_pdg").c_str(), "PDG", w1 - w0, w0, w1, d1 - d0, d0, d1);
          }

          for (size_t w = w0; w < w1; ++w) {
            auto const& raw = fTrainingDataAlg.wireData(w);
            for (size_t d = d0; d < d1; ++d) {
              rawHist->Fill(w, d, raw[d]);
            }

            if (saveSim) {
              auto const& edep = fTrainingDataAlg.wireEdep(w);
              for (size_t d = d0; d < d1; ++d) {
                depHist->Fill(w, d, edep[d]);
              }

              auto const& pdg = fTrainingDataAlg.wirePdg(w);
              for (size_t d = d0; d < d1; ++d) {
                pdgHist->Fill(w, d, pdg[d]);
              }
            }
          }

          writeAndDelete(rawHist);
          writeAndDelete(depHist);
          writeAndDelete(pdgHist);
        }
        else if (fDumpToNumpy) {
          for (size_t w = w0; w < w1; ++w) {
            int w_start = w - fPatch_size_w / 2;
            int w_stop = w_start + fPatch_size_w;
            if (w_start < int(w0) || w_start > int(w1)) continue;
            if (w_stop < int(w0) || w_stop > int(w1)) continue;
            auto const& pdg = fTrainingDataAlg.wirePdg(w);
            auto const& deposit = fTrainingDataAlg.wireEdep(w);
            auto const& raw = fTrainingDataAlg.wireData(w);
            for (size_t d = d0; d < d1; ++d) {
              int d_start = d - fPatch_size_d / 2;
              int d_stop = d_start + fPatch_size_d;
              if (d_start < int(d0) || d_start > int(d1)) continue;
              if (d_stop < int(d0) || d_stop > int(d1)) continue;
              int y0 = 0, y1 = 0, y2 = 0, y3 = 0;
              if (deposit[d] < 2e-5 || raw[d] < 0.05) { //empty pixel
                y3 = 1;
                ++nNone;
                if (flat.fire() > fNone) continue;
                ++nNone_sel;
              }
              else if ((pdg[d] & 0x0FFF) == 11) { //shower
                y1 = 1;
                ++nEm;
                if ((pdg[d] & 0xF000) == 0x2000) { //Michel
                  y2 = 1;
                  ++nMichel;
                  if (flat.fire() > fMichel) continue;
                  ++nMichel_sel;
                }
                else {
                  if (flat.fire() > fEm) continue;
                  ++nEm_sel;
                }
              }
              else { //track
                y0 = 1;
                ++nTrk;
                //Check if the track appears to be stopping
                int nPxlw0 = 0, nPxlw1 = 0, nPxld0 = 0, nPxld1 = 0;
                //Try to fit a track
                std::vector<double> wfit, dfit, chgfit;
                double total_trkchg = 0;
                for (int ww = w_start; ww < w_stop; ++ww) {
                  auto const& pdg1 = fTrainingDataAlg.wirePdg(ww);
                  auto const& deposit1 = fTrainingDataAlg.wireEdep(ww);
                  auto const& raw1 = fTrainingDataAlg.wireData(ww);
                  for (int dd = d_start; dd < d_stop; ++dd) {
                    if (deposit1[dd] < 2e-5 || raw1[dd] < 0.05) continue; //empty pixel
                    if ((pdg1[dd] & 0x0FFF) != 11) {                      //track pixel
                      wfit.push_back(ww);
                      dfit.push_back(dd);
                      chgfit.push_back(raw1[dd]);
                      total_trkchg += raw1[dd];
                      if (ww - w_start < 3) ++nPxlw0;
                      if (w_stop - ww - 1 < 3) ++nPxlw1;
                      if (dd - d_start < 3) ++nPxld0;
                      if (d_stop - dd - 1 < 3) ++nPxld1;
                    }
                  }
                }
                //Fit track pixels
                double fit_trkchg = 0;
                double parm[2] = {};
                if (!wfit.empty()) {
                  if (!WeightedFit(wfit.size(), wfit, dfit, chgfit, &parm[0])) {
                    for (size_t j = 0; j < wfit.size(); ++j) {
                      if (std::abs((dfit[j] - (parm[0] + wfit[j] * parm[1])) * cos(atan(parm[1]))) <
                          3) {
                        fit_trkchg += chgfit[j];
                      }
                    }
                  }
                  else if (!WeightedFit(dfit.size(), dfit, wfit, chgfit, &parm[0])) {
                    for (size_t j = 0; j < dfit.size(); ++j) {
                      if (std::abs((wfit[j] - (parm[0] + dfit[j] * parm[1])) * cos(atan(parm[1]))) <
                          3) {
                        fit_trkchg += chgfit[j];
                      }
                    }
                  }
                }
                if (((nPxlw0) && (!nPxlw1) && (!nPxld0) && (!nPxld1)) ||
                    ((!nPxlw0) && (nPxlw1) && (!nPxld0) && (!nPxld1)) ||
                    ((!nPxlw0) && (!nPxlw1) && (nPxld0) && (!nPxld1)) ||
                    ((!nPxlw0) && (!nPxlw1) && (!nPxld0) && (nPxld1))) { //stopping track
                  if (flat.fire() > fStopTrk) continue;
                  ++nStopTrk_sel;
                }
                else if (total_trkchg && fit_trkchg / total_trkchg > 0.9) { //clean track
                  if (flat.fire() > fCleanTrk) continue;
                  ++nCleanTrk_sel;
                }
                else {
                  if (flat.fire() > fTrk) continue;
                  ++nTrk_sel;
                }
              }
              c2numpy_uint32(&npywriter, event.id().run());
              c2numpy_uint32(&npywriter, event.id().subRun());
              c2numpy_uint32(&npywriter, event.id().event());
              c2numpy_uint8(&npywriter, fSelectedTPC[i]);
              c2numpy_uint8(&npywriter, fSelectedPlane[v]);
              c2numpy_uint16(&npywriter, w);
              c2numpy_uint16(&npywriter, d);
              c2numpy_uint8(&npywriter, y0);
              c2numpy_uint8(&npywriter, y1);
              c2numpy_uint8(&npywriter, y2);
              c2numpy_uint8(&npywriter, y3);
              for (int ww = w_start; ww < w_stop; ++ww) {
                auto const& raw1 = fTrainingDataAlg.wireData(ww);
                for (int dd = d_start; dd < d_stop; ++dd) {
                  c2numpy_float32(&npywriter, raw1[dd]);
                }
              }
            }
          }
        }
        else {
          std::ostringstream ss1;
          ss1 << fOutTextFilePath << "/raw_" << os.str() << "_tpc_" << fSelectedTPC[i] << "_view_"
              << fSelectedPlane[v];

          std::ofstream fout_raw, fout_deposit, fout_pdg;

          fout_raw.open(ss1.str() + ".raw");
          if (saveSim) {
            fout_deposit.open(ss1.str() + ".deposit");
            fout_pdg.open(ss1.str() + ".pdg");
          }

          for (size_t w = w0; w < w1; ++w) {
            auto const& raw = fTrainingDataAlg.wireData(w);
            for (size_t d = d0; d < d1; ++d) {
              fout_raw << raw[d] << " ";
            }
            fout_raw << std::endl;

            if (saveSim) {
              auto const& edep = fTrainingDataAlg.wireEdep(w);
              for (size_t d = d0; d < d1; ++d) {
                fout_deposit << edep[d] << " ";
              }
              fout_deposit << std::endl;

              auto const& pdg = fTrainingDataAlg.wirePdg(w);
              for (size_t d = d0; d < d1; ++d) {
                fout_pdg << pdg[d] << " ";
              }
              fout_pdg << std::endl;
            }
          }

          fout_raw.close();
          if (saveSim) {
            fout_deposit.close();
            fout_pdg.close();
          }
        }
      }

  } // PointIdTrainingData::analyze()

  //-----------------------------------------------------------------------
  int PointIdTrainingData::WeightedFit(int n,
                                       std::vector<double> const& x,
                                       std::vector<double> const& y,
                                       std::vector<double> const& w,
                                       double* params) const
  {
    Double_t sumx = 0.;
    Double_t sumx2 = 0.;
    Double_t sumy = 0.;
    Double_t sumy2 = 0.;
    Double_t sumxy = 0.;
    Double_t sumw = 0.;
    Double_t eparams[2] = {};

    for (Int_t i = 0; i < n; i++) {
      sumx += x[i] * w[i];
      sumx2 += x[i] * x[i] * w[i];
      sumy += y[i] * w[i];
      sumy2 += y[i] * y[i] * w[i];
      sumxy += x[i] * y[i] * w[i];
      sumw += w[i];
    }

    if (sumx2 * sumw - sumx * sumx == 0.) return 1;
    if (sumx2 - sumx * sumx / sumw == 0.) return 1;

    params[0] = (sumy * sumx2 - sumx * sumxy) / (sumx2 * sumw - sumx * sumx);
    params[1] = (sumxy - sumx * sumy / sumw) / (sumx2 - sumx * sumx / sumw);

    eparams[0] = sumx2 * (sumx2 * sumw - sumx * sumx);
    eparams[1] = (sumx2 - sumx * sumx / sumw);

    if (eparams[0] < 0. || eparams[1] < 0.) return 1;

    eparams[0] = sqrt(eparams[0]) / (sumx2 * sumw - sumx * sumx);
    eparams[1] = sqrt(eparams[1]) / (sumx2 - sumx * sumx / sumw);

    return 0;
  }

  DEFINE_ART_MODULE(PointIdTrainingData)

}
