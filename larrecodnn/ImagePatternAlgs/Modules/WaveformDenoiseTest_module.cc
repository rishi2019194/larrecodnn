//
// This is for testing the 1d denoising AE
//
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IWaveformDenoise.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IWaveformRecog.h"

#include "c2numpy.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;

struct CnnROI {
  unsigned int start;
  unsigned int end;
};

namespace nnet {
  class WaveformDenoiseTest;
}

class nnet::WaveformDenoiseTest : public art::EDAnalyzer {

public:
  explicit WaveformDenoiseTest(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  WaveformDenoiseTest(WaveformDenoiseTest const&) = delete;
  WaveformDenoiseTest(WaveformDenoiseTest&&) = delete;
  WaveformDenoiseTest& operator=(WaveformDenoiseTest const&) = delete;
  WaveformDenoiseTest& operator=(WaveformDenoiseTest&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  //void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;
  void endJob() override;

private:
  std::string fDumpWaveformsFileName;
  std::string fDumpCleanSignalFileName;
  std::string fDumpDenoisedFileName;
  unsigned int fDisplayWindowSize;

  std::string fSimChannelLabel;             ///< module that made simchannels
  std::string fDigitModuleLabel;            ///< module that made digits
  std::string fCleanSignalDigitModuleLabel; ///< module that made the signal-only digits

  std::string fPlaneToDump;
  std::string fCollectionPlaneLabel;
  art::ServiceHandle<geo::Geometry> fgeom;

  std::vector<std::unique_ptr<wavrec_tool::IWaveformRecog>> fWaveformRecogToolVec;
  std::vector<std::unique_ptr<wavdenoise_tool::IWaveformDenoise>> fWaveformDenoiseToolVec;
  int fNPlanes;

  std::ofstream _fileout1;

  c2numpy_writer npywriter;
  c2numpy_writer npywriter2;
  c2numpy_writer npywriter3;
};

//-----------------------------------------------------------------------
nnet::WaveformDenoiseTest::WaveformDenoiseTest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fDumpWaveformsFileName(p.get<std::string>("DumpWaveformsFileName", "dumpwaveforms"))
  , fDumpCleanSignalFileName(p.get<std::string>("CleanSignalFileName", "dumpcleansignal"))
  , fDumpDenoisedFileName(p.get<std::string>("DenoisedFileName", "dumpdenoised"))
  , fDisplayWindowSize(p.get<unsigned int>("DisplayWindowSize", 500))
  , fSimChannelLabel(p.get<std::string>("SimChannelLabel", "elecDrift"))
  , fDigitModuleLabel(p.get<std::string>("DigitModuleLabel", "simWire"))
  , fCleanSignalDigitModuleLabel(
      p.get<std::string>("CleanSignalDigitModuleLabel", "simWire:signal"))
  , fPlaneToDump(p.get<std::string>("PlaneToDump", "U"))
  , fCollectionPlaneLabel(p.get<std::string>("CollectionPlaneLabel", "Z"))
{

  if (fDigitModuleLabel.empty() && fCleanSignalDigitModuleLabel.empty()) {
    throw cet::exception("WaveformDenoiseTest")
      << "Both DigitModuleLabel and CleanSignalModuleLabel are empty";
  }

  fNPlanes = fgeom->Nplanes();

  // Signal/Noise waveform recognition tool
  fWaveformRecogToolVec.reserve(fNPlanes);
  auto const tool_psets1 = p.get<std::vector<fhicl::ParameterSet>>("WaveformRecogs");
  for (auto const& pset : tool_psets1) {
    fWaveformRecogToolVec.push_back(art::make_tool<wavrec_tool::IWaveformRecog>(pset));
  }

  // AE based waveform denoising tool
  fWaveformDenoiseToolVec.reserve(fNPlanes);
  auto const tool_psets2 = p.get<std::vector<fhicl::ParameterSet>>("WaveformDenoisers");
  for (auto const& pset : tool_psets2) {
    fWaveformDenoiseToolVec.push_back(art::make_tool<wavdenoise_tool::IWaveformDenoise>(pset));
  }
}

//-----------------------------------------------------------------------
void nnet::WaveformDenoiseTest::beginJob()
{

  c2numpy_init(&npywriter, fDumpWaveformsFileName, 50000);
  c2numpy_addcolumn(&npywriter, "evt", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "chan", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "roi", C2NUMPY_UINT16);
  c2numpy_addcolumn(&npywriter, "view", (c2numpy_type)((int)C2NUMPY_STRING + 1));
  for (unsigned int i = 0; i < fDisplayWindowSize; i++) {
    std::ostringstream name;
    name << "tck_" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT16);
  }

  // ... this is for storing the clean signal (no noise) waveform
  c2numpy_init(&npywriter2, fDumpCleanSignalFileName, 50000);

  for (unsigned int i = 0; i < fDisplayWindowSize; i++) {
    std::ostringstream name;
    name << "tck_" << i;
    c2numpy_addcolumn(&npywriter2, name.str().c_str(), C2NUMPY_INT16);
  }

  // ... this is for storing the denoised raw waveform
  c2numpy_init(&npywriter3, fDumpDenoisedFileName, 50000);

  for (unsigned int i = 0; i < fDisplayWindowSize; i++) {
    std::ostringstream name;
    name << "tck_" << i;
    c2numpy_addcolumn(&npywriter3, name.str().c_str(), C2NUMPY_FLOAT32);
  }
}

//-----------------------------------------------------------------------
void nnet::WaveformDenoiseTest::endJob()
{
  c2numpy_close(&npywriter);
  c2numpy_close(&npywriter2);
  c2numpy_close(&npywriter3);
}

//-----------------------------------------------------------------------
void nnet::WaveformDenoiseTest::analyze(art::Event const& evt)
{
  cout << "Event "
       << " " << evt.id().run() << " " << evt.id().subRun() << " " << evt.id().event() << endl;

  // ... Read in the digit List object(s).
  auto digitVecHandle = evt.getValidHandle<std::vector<raw::RawDigit>>(fDigitModuleLabel);

  // ... Read in the signal-only digit List object(s).
  auto digitVecHandle2 =
    evt.getValidHandle<std::vector<raw::RawDigit>>(fCleanSignalDigitModuleLabel);

  if (!digitVecHandle->size() || !digitVecHandle2->size()) {
    throw cet::exception("NoiseWaveformRoiAnalyzer")
      << "At least one of the raw digits lists is empty.";
  }

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  // ... Use the handle to get a particular (0th) element of collection.
  std::vector<raw::RawDigit> const& rawDigits = *digitVecHandle;
  raw::RawDigit const& digitVec0 = rawDigits[0];
  unsigned int dataSize = digitVec0.Samples(); //size of raw data vectors
  if (dataSize != detProp.ReadOutWindowSize()) {
    throw cet::exception("NoiseWaveformRoiAnalyzer") << "Bad dataSize: " << dataSize;
  }
  std::vector<raw::RawDigit> const& rawDigits2 = *digitVecHandle2;
  raw::RawDigit const& digitVec20 = rawDigits2[0];
  unsigned int dataSize2 = digitVec20.Samples(); //size of raw data vectors
  if (dataSize != dataSize2) {
    throw cet::exception("NoiseWaveformRoiAnalyzer")
      << "RawDigits from the 2 data products have different dataSizes: " << dataSize << "not eq to"
      << dataSize2;
  }

  // ... Build a map from channel number -> rawdigitVec
  std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> rawdigitMap;
  raw::ChannelID_t chnum = raw::InvalidChannelID; // channel number
  for (size_t rdIter = 0; rdIter < rawDigits.size(); ++rdIter) {
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
    chnum = digitVec->Channel();
    if (chnum == raw::InvalidChannelID) continue;
    if (geo::PlaneGeo::ViewName(fgeom->View(chnum)) != fPlaneToDump[0]) continue;
    rawdigitMap[chnum] = digitVec;
  }
  std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> rawdigitMap2;
  raw::ChannelID_t chnum2 = raw::InvalidChannelID; // channel number
  for (size_t rdIter = 0; rdIter < rawDigits2.size(); ++rdIter) {
    art::Ptr<raw::RawDigit> digitVec2(digitVecHandle2, rdIter);
    chnum2 = digitVec2->Channel();
    if (chnum2 == raw::InvalidChannelID) continue;
    if (geo::PlaneGeo::ViewName(fgeom->View(chnum2)) != fPlaneToDump[0]) continue;
    rawdigitMap2[chnum2] = digitVec2;
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ... First find all channels that have non-zero raw digits and erase them
  //     from the raw digit maps
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  std::vector<raw::ChannelID_t> signalchannels;

  for (std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>>::iterator iter = rawdigitMap2.begin();
       iter != rawdigitMap2.end();
       ++iter) {
    std::vector<short> rawadc(dataSize);  // vector to hold uncompressed adc values later
    std::vector<short> adcvec2(dataSize); // vector to pedestal-subtracted adc values

    auto rawdig2 = iter->second;
    raw::Uncompress(rawdig2->ADCs(), rawadc, rawdig2->GetPedestal(), rawdig2->Compression());
    for (size_t j = 0; j < rawadc.size(); ++j) {
      adcvec2[j] = rawadc[j] - rawdig2->GetPedestal();
    }

    auto itnz = find_if(adcvec2.begin(), adcvec2.end(), [](auto x) { return x != 0; });
    if (itnz != adcvec2.end()) signalchannels.push_back(iter->first);
  }

  for (size_t ich = 0; ich < signalchannels.size(); ++ich) {

    // .. get signal-containing channel number
    auto ch1 = signalchannels[ich];
    if (ch1 == raw::InvalidChannelID) continue;
    if (geo::PlaneGeo::ViewName(fgeom->View(ch1)) != fPlaneToDump[0]) continue;

    std::vector<short> rawadc(dataSize);  // vector to hold uncompressed adc values later
    std::vector<short> adcvec(dataSize);  // vector to hold zero-padded full waveform
    std::vector<short> adcvec2(dataSize); // vector to hold zero-padded full signal-only waveform
    std::vector<float> inputsignal(dataSize); // float version of adcvec

    auto search = rawdigitMap.find(ch1);
    if (search == rawdigitMap.end()) continue;
    auto rawdig = (*search).second;
    raw::Uncompress(rawdig->ADCs(), rawadc, rawdig->GetPedestal(), rawdig->Compression());
    for (size_t j = 0; j < rawadc.size(); ++j) {
      adcvec[j] = rawadc[j] - rawdig->GetPedestal();
      inputsignal[j] = adcvec[j];
    }

    auto search2 = rawdigitMap2.find(ch1);
    if (search2 == rawdigitMap2.end()) continue;
    auto rawdig2 = (*search2).second;
    raw::Uncompress(rawdig2->ADCs(), rawadc, rawdig2->GetPedestal(), rawdig2->Compression());
    for (size_t i = 0; i < rawadc.size(); ++i) {
      adcvec2[i] = rawadc[i] - rawdig2->GetPedestal();
    }

    // ... use waveform recognition CNN to perform inference on each window
    std::vector<bool> inroi(dataSize, false);
    inroi = fWaveformRecogToolVec[fgeom->View(ch1)]->findROI(inputsignal);

    auto itnf = find_if(inroi.begin(), inroi.end(), [](auto x) { return x; });
    if (itnf == inroi.end()) continue;

    CnnROI roi;
    std::vector<CnnROI> cnn_rois;

    bool is_roi = false;
    bool was_roi = false;

    for (size_t itck = 0; itck < inroi.size() - 1; ++itck) {
      is_roi = inroi[itck];
      if (is_roi && !was_roi) { roi.start = itck; }
      else if (!is_roi && was_roi) {
        roi.end = itck - 1;
        cnn_rois.push_back(roi);
      }
      was_roi = is_roi;
    }

    for (size_t i = 0; i < cnn_rois.size(); ++i) {
      std::vector<short> wavraw(fDisplayWindowSize, 0);
      std::vector<short> wavcln(fDisplayWindowSize, 0);
      std::vector<float> wavdns(fDisplayWindowSize, 0);

      unsigned int tcka = cnn_rois[i].start;
      unsigned int tckb = cnn_rois[i].end;
      unsigned int ntcks_roi = tckb - tcka + 1;
      unsigned int ntcks_dis = ntcks_roi;
      if (ntcks_dis > fDisplayWindowSize) {
        tckb = tcka + fDisplayWindowSize - 1;
        ntcks_dis = tckb - tcka + 1;
      }

      unsigned jtck = 0;
      for (unsigned int itck = tcka; itck < tckb; ++itck) {
        wavraw[jtck] = adcvec[itck];
        wavcln[jtck] = adcvec2[itck];
        jtck++;
      }

      std::vector<float> wavinp(ntcks_roi, 0);
      std::vector<float> wavout(ntcks_roi, 0);

      jtck = 0;
      for (unsigned int itck = tcka; itck < tcka + ntcks_roi; ++itck) {
        wavinp[jtck] = float(adcvec[itck]);
        jtck++;
      }
      wavout = fWaveformDenoiseToolVec[fgeom->View(ch1)]->denoiseWaveform(wavinp);

      for (unsigned int itck = 0; itck < ntcks_dis; ++itck) {
        wavdns[itck] = wavout[itck];
      }

      c2numpy_uint32(&npywriter, evt.id().event());
      c2numpy_uint32(&npywriter, ch1);
      c2numpy_uint16(&npywriter, i); // roi
      c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(ch1)).c_str());

      for (unsigned int j = 0; j < fDisplayWindowSize; ++j) {
        c2numpy_int16(&npywriter, wavraw[j]);
      }
      for (unsigned int j = 0; j < fDisplayWindowSize; ++j) {
        c2numpy_int16(&npywriter2, wavcln[j]);
      }
      for (unsigned int j = 0; j < fDisplayWindowSize; ++j) {
        c2numpy_float32(&npywriter3, wavdns[j]);
      }
    }
  }
}
DEFINE_ART_MODULE(nnet::WaveformDenoiseTest)
