//////////////////////////////////////////////
//
//  Module to dump noise waveforms for 1DCNN
//  training
//
//  mwang@fnal.gov
//
//////////////////////////////////////////////

#include <random>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
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
#include "lardataobj/Simulation/SimChannel.h"

#include "CLHEP/Random/RandFlat.h"
#include "c2numpy.h"
#include "nurandom/RandomUtils/NuRandomService.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;

namespace nnet {
  class NoiseWaveformDump;
}

class nnet::NoiseWaveformDump : public art::EDAnalyzer {

public:
  explicit NoiseWaveformDump(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  NoiseWaveformDump(NoiseWaveformDump const&) = delete;
  NoiseWaveformDump(NoiseWaveformDump&&) = delete;
  NoiseWaveformDump& operator=(NoiseWaveformDump const&) = delete;
  NoiseWaveformDump& operator=(NoiseWaveformDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  //void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;
  void endJob() override;

private:
  std::string fDumpWaveformsFileName;

  std::string fSimulationProducerLabel; ///< producer that tracked simulated part. through detector
  std::string fSimChannelLabel;         ///< module that made simchannels
  std::string fDigitModuleLabel;        ///< module that made digits
  bool fUseFullWaveform;
  unsigned int fShortWaveformSize;

  std::string fPlaneToDump;
  int fMaxNoiseChannelsPerEvent;
  art::ServiceHandle<geo::Geometry> fgeom;

  CLHEP::RandFlat fRandFlat;

  c2numpy_writer npywriter;
};

// Create the random number generator
namespace {
  std::string const instanceName = "NoiseWaveformDump";
}

//-----------------------------------------------------------------------
nnet::NoiseWaveformDump::NoiseWaveformDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fDumpWaveformsFileName(p.get<std::string>("DumpWaveformsFileName", "dumpwaveforms"))
  , fSimulationProducerLabel(p.get<std::string>("SimulationProducerLabel", "larg4Main"))
  , fSimChannelLabel(p.get<std::string>("SimChannelLabel", "elecDrift"))
  , fDigitModuleLabel(p.get<std::string>("DigitModuleLabel", "simWire"))
  , fUseFullWaveform(p.get<bool>("UseFullWaveform", true))
  , fShortWaveformSize(p.get<unsigned int>("ShortWaveformSize"))
  , fPlaneToDump(p.get<std::string>("PlaneToDump"))
  , fMaxNoiseChannelsPerEvent(p.get<int>("MaxNoiseChannelsPerEvent"))
  , fRandFlat{createEngine(
      art::ServiceHandle<rndm::NuRandomService> {}
      -> declareEngine(instanceName, p, "SeedForNoiseWaveformDump"),
      "HepJamesRandom",
      instanceName)}
{
  if (std::getenv("CLUSTER") && std::getenv("PROCESS")) {
    fDumpWaveformsFileName +=
      string(std::getenv("CLUSTER")) + "-" + string(std::getenv("PROCESS")) + "-";
  }

  if (fDigitModuleLabel.empty()) {
    throw cet::exception("NoiseWaveformDump") << "DigitModuleLabel is empty";
  }
}

//-----------------------------------------------------------------------
void nnet::NoiseWaveformDump::beginJob()
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();

  c2numpy_init(&npywriter, fDumpWaveformsFileName, 50000);
  c2numpy_addcolumn(&npywriter, "evt", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "chan", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "view", (c2numpy_type)((int)C2NUMPY_STRING + 1));
  c2numpy_addcolumn(&npywriter, "ntrk", C2NUMPY_UINT16);

  for (unsigned int i = 0; i < 5; i++) {
    std::ostringstream name;

    name.str("");
    name << "tid" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT32);

    name.str("");
    name << "pdg" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT32);

    name.str("");
    name << "gen" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), (c2numpy_type)((int)C2NUMPY_STRING + 6));

    name.str("");
    name << "pid" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), (c2numpy_type)((int)C2NUMPY_STRING + 7));

    name.str("");
    name << "edp" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_FLOAT32);

    name.str("");
    name << "nel" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_UINT32);

    name.str("");
    name << "sti" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_UINT16);

    name.str("");
    name << "stf" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_UINT16);
  }

  for (unsigned int i = 0;
       i < (fUseFullWaveform ? detProp.ReadOutWindowSize() : fShortWaveformSize);
       i++) {
    std::ostringstream name;
    name << "tck_" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT16);
  }
}

//-----------------------------------------------------------------------
void nnet::NoiseWaveformDump::endJob()
{
  c2numpy_close(&npywriter);
}

//-----------------------------------------------------------------------
void nnet::NoiseWaveformDump::analyze(art::Event const& evt)
{
  cout << "Event " << evt.id().run() << " " << evt.id().subRun() << " " << evt.id().event() << endl;

  // ... Read in the digit List object(s).
  auto digitVecHandle = evt.getValidHandle<std::vector<raw::RawDigit>>(fDigitModuleLabel);

  std::cout << " !!!!! Size of digitVecHandle: " << digitVecHandle->size() << std::endl;
  if (!digitVecHandle->size()) return;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  // ... Use the handle to get a particular (0th) element of collection.
  std::vector<raw::RawDigit> const& rawDigits = *digitVecHandle;
  raw::RawDigit const& digitVec0 = rawDigits[0];
  unsigned int dataSize = digitVec0.Samples(); //size of raw data vectors
  if (dataSize != detProp.ReadOutWindowSize()) {
    throw cet::exception("NoiseWaveformDumpNoiseWaveformDump") << "Bad dataSize: " << dataSize;
  }

  // ... Build a map from channel number -> rawdigitVec
  std::map<raw::ChannelID_t, raw::RawDigit const*> rawdigitMap;
  raw::ChannelID_t chnum = raw::InvalidChannelID; // channel number
  for (size_t rdIter = 0; rdIter < rawDigits.size(); ++rdIter) {
    chnum = rawDigits[rdIter].Channel();
    if (chnum == raw::InvalidChannelID) continue;
    if (geo::PlaneGeo::ViewName(fgeom->View(chnum)) != fPlaneToDump[0]) continue;
    rawdigitMap[chnum] = &rawDigits[rdIter];
  }

  // ... Read in sim channel list
  auto simChannelHandle = evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);
  std::cout << " !!!!! Size of simChannelHandle: " << simChannelHandle->size() << std::endl;

  if (simChannelHandle->size() > 0) {

    // ... Loop over simChannels to find all signal channels and erase them
    //     from the rawdigitMap
    for (auto const& channel : (*simChannelHandle)) {

      // .. get simChannel channel number
      const raw::ChannelID_t ch1 = channel.Channel();
      if (ch1 == raw::InvalidChannelID) continue;
      if (geo::PlaneGeo::ViewName(fgeom->View(ch1)) != fPlaneToDump[0]) continue;

      bool hasEnergyDeposit = false;

      // ... Loop over all ticks with ionization energy deposited
      auto const& timeSlices = channel.TDCIDEMap();
      for (auto const& [tpctime, energyDeposits] : timeSlices) {

        unsigned int tdctick = static_cast<unsigned int>(clockData.TPCTDC2Tick(double(tpctime)));
        if (tdctick < 0 || tdctick > (dataSize - 1)) continue;

        // ... Loop over all energy depositions in this tick
        for (auto const& energyDeposit : energyDeposits) {
          if ((energyDeposit.energy > 0) || (energyDeposit.numElectrons > 0)) {
            hasEnergyDeposit = true;
            break;
          }
        }
        if (hasEnergyDeposit) break;
      }
      if (hasEnergyDeposit) rawdigitMap.erase(ch1);

    } // loop over SimChannels
  }

  // .. Now construct noise channel vector from updated rawdigitMap
  std::vector<raw::ChannelID_t> noisechannels;
  for (std::map<raw::ChannelID_t, raw::RawDigit const*>::iterator iter = rawdigitMap.begin();
       iter != rawdigitMap.end();
       ++iter) {
    noisechannels.push_back(iter->first);
  }

  std::cout << " !!!!! size of noisechannels: " << noisechannels.size() << std::endl;
  if (!noisechannels.size()) return;

  //save noise
  int noisechancount = 0;

  // .. create a vector for shuffling the wire channel indices
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::vector<size_t> randigitmap;
  for (size_t i = 0; i < noisechannels.size(); ++i)
    randigitmap.push_back(i);
  std::shuffle(randigitmap.begin(), randigitmap.end(), std::mt19937(seed));

  std::string dummystr6 = "none  ";
  std::string dummystr7 = "none   ";

  for (size_t rdIter = 0; rdIter < noisechannels.size(); ++rdIter) {

    if (noisechancount == fMaxNoiseChannelsPerEvent) break;

    size_t ranIdx = randigitmap[rdIter];
    auto digitVec = rawdigitMap[noisechannels[ranIdx]];

    std::vector<short> rawadc(dataSize); // vector to hold uncompressed adc values later
    std::vector<short> adcvec(dataSize); // vector to pedestal-subtracted adc values
    raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
    for (size_t j = 0; j < rawadc.size(); ++j) {
      adcvec[j] = rawadc[j] - digitVec->GetPedestal();
    }
    c2numpy_uint32(&npywriter, evt.id().event());
    c2numpy_uint32(&npywriter, digitVec->Channel());
    c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(digitVec->Channel())).c_str());

    c2numpy_uint16(&npywriter, 0); //number of peaks
    for (unsigned int i = 0; i < 5; ++i) {
      c2numpy_int32(&npywriter, 0);
      c2numpy_int32(&npywriter, 0);
      c2numpy_string(&npywriter, dummystr6.c_str());
      c2numpy_string(&npywriter, dummystr7.c_str());
      c2numpy_float32(&npywriter, 0.);
      c2numpy_uint32(&npywriter, 0);
      c2numpy_uint16(&npywriter, 0);
      c2numpy_uint16(&npywriter, 0);
    }

    if (fUseFullWaveform) {
      for (unsigned int itck = 0; itck < dataSize; ++itck) {
        c2numpy_int16(&npywriter, adcvec[itck]);
      }
    }
    else {
      int start_tick = int((dataSize - fShortWaveformSize) * fRandFlat.fire(0, 1));
      for (unsigned int itck = start_tick; itck < (start_tick + fShortWaveformSize); ++itck) {
        c2numpy_int16(&npywriter, adcvec[itck]);
      }
    }

    ++noisechancount;
  }
  std::cout << "Total number of noise channels " << noisechancount << std::endl;
}
DEFINE_ART_MODULE(nnet::NoiseWaveformDump)
