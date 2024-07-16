//////////////////////////////////////////////
//
//  Module to dump raw data for ROI training
//
//  mwang@fnal.gov
//  tjyang@fnal.gov
//  wwu@fnal.gov
//
//  This version:
//    a) uses RawDigit or WireProducer
//    b) saves full waveform or short waveform
//       - for short waveform, picks signal
//         associated with largest energy
//         deposit & saves only tdcmin/tdcmaxes
//         that are within window extents
//
//////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "CLHEP/Random/RandFlat.h"
#include "c2numpy.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <map>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::string;

struct WireSigInfo {
  int pdgcode;
  std::string genlab;
  std::string procid;
  unsigned int tdcmin;
  unsigned int tdcmax;
  int numel;
  double edep;
};

namespace nnet {
  class RawWaveformDump;
}

class nnet::RawWaveformDump : public art::EDAnalyzer {

public:
  explicit RawWaveformDump(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  RawWaveformDump(RawWaveformDump const&) = delete;
  RawWaveformDump(RawWaveformDump&&) = delete;
  RawWaveformDump& operator=(RawWaveformDump const&) = delete;
  RawWaveformDump& operator=(RawWaveformDump&&) = delete;

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
  std::string fWireProducerLabel;
  bool fUseFullWaveform;
  unsigned int fShortWaveformSize;

  std::string fSelectGenLabel;
  std::string fSelectProcID;
  int fSelectPDGCode;
  std::string fPlaneToDump;
  double fMinParticleEnergyGeV;
  double fMinEnergyDepositedMeV;
  int fMinNumberOfElectrons;
  int fMaxNumberOfElectrons;
  bool fSaveSignal;
  int fMaxNoiseChannelsPerEvent;
  std::string fCollectionPlaneLabel;
  art::ServiceHandle<geo::Geometry> fgeom;
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;

  CLHEP::RandFlat fRandFlat;

  c2numpy_writer npywriter;
};

//-----------------------------------------------------------------------
struct genFinder {
private:
  typedef std::pair<int, std::string> track_id_to_string;
  std::vector<track_id_to_string> track_id_map;
  std::set<std::string> generator_names;
  bool isSorted = false;

public:
  void sort_now()
  {
    std::sort(this->track_id_map.begin(),
              this->track_id_map.end(),
              [](const auto& a, const auto& b) { return (a.first < b.first); });
    isSorted = true;
  }
  void add(const int& track_id, const std::string& gname)
  {
    this->track_id_map.push_back(std::make_pair(track_id, gname));
    generator_names.emplace(gname);
    isSorted = false;
  }
  bool has_gen(std::string gname) { return static_cast<bool>(generator_names.count(gname)); };
  std::string get_gen(int tid)
  {
    if (!isSorted) { this->sort_now(); }
    return std::lower_bound(track_id_map.begin(),
                            track_id_map.end(),
                            tid,
                            [](const auto& a, const auto& b) { return (a.first < b); })
      ->second;
  };
};

// Create the random number generator
namespace {
  std::string const instanceName = "RawWaveformDump";
}

//-----------------------------------------------------------------------
nnet::RawWaveformDump::RawWaveformDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fDumpWaveformsFileName(p.get<std::string>("DumpWaveformsFileName", "dumpwaveforms"))
  , fSimulationProducerLabel(p.get<std::string>("SimulationProducerLabel", "larg4Main"))
  , fSimChannelLabel(p.get<std::string>("SimChannelLabel", "elecDrift"))
  , fDigitModuleLabel(p.get<std::string>("DigitModuleLabel", "simWire"))
  , fWireProducerLabel(p.get<std::string>("WireProducerLabel"))
  , fUseFullWaveform(p.get<bool>("UseFullWaveform", true))
  , fShortWaveformSize(p.get<unsigned int>("ShortWaveformSize"))
  , fSelectGenLabel(p.get<std::string>("SelectGenLabel", "ANY"))
  , fSelectProcID(p.get<std::string>("SelectProcID", "ANY"))
  , fSelectPDGCode(p.get<int>("SelectPDGCode", 0))
  , fPlaneToDump(p.get<std::string>("PlaneToDump"))
  , fMinParticleEnergyGeV(p.get<double>("MinParticleEnergyGeV", 0.))
  , fMinEnergyDepositedMeV(p.get<double>("MinEnergyDepositedMeV", 0.))
  , fMinNumberOfElectrons(p.get<int>("MinNumberOfElectrons", 600))
  , fMaxNumberOfElectrons(p.get<int>("MaxNumberOfElectrons", -1))
  , fSaveSignal(p.get<bool>("SaveSignal", true))
  , fMaxNoiseChannelsPerEvent(p.get<int>("MaxNoiseChannelsPerEvent"))
  , fCollectionPlaneLabel(p.get<std::string>("CollectionPlaneLabel"))
  , fRandFlat{createEngine(
      art::ServiceHandle<rndm::NuRandomService> {}
      -> declareEngine(instanceName, p, "SeedForRawWaveformDump"),
      "HepJamesRandom",
      instanceName)}
{
  if (std::getenv("CLUSTER") && std::getenv("PROCESS")) {
    fDumpWaveformsFileName +=
      string(std::getenv("CLUSTER")) + "-" + string(std::getenv("PROCESS")) + "-";
  }

  if (fDigitModuleLabel.empty() && fWireProducerLabel.empty()) {
    throw cet::exception("RawWaveformDump")
      << "Both DigitModuleLabel and WireProducerLabel are empty";
  }

  if ((!fDigitModuleLabel.empty()) && (!fWireProducerLabel.empty())) {
    throw cet::exception("RawWaveformDump")
      << "Only one of DigitModuleLabel and WireProducerLabel should be set";
  }
}

//-----------------------------------------------------------------------
void nnet::RawWaveformDump::beginJob()
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
void nnet::RawWaveformDump::endJob()
{
  c2numpy_close(&npywriter);
}

//-----------------------------------------------------------------------
void nnet::RawWaveformDump::analyze(art::Event const& evt)
{
  cout << "Event "
       << " " << evt.id().run() << " " << evt.id().subRun() << " " << evt.id().event() << endl;

  std::unique_ptr<genFinder> gf(new genFinder());

  // ... Read in the digit List object(s).
  art::Handle<std::vector<raw::RawDigit>> digitVecHandle;
  std::vector<art::Ptr<raw::RawDigit>> rawdigitlist;
  if (evt.getByLabel(fDigitModuleLabel, digitVecHandle)) {
    art::fill_ptr_vector(rawdigitlist, digitVecHandle);
  }

  // ... Read in the wire List object(s).
  art::Handle<std::vector<recob::Wire>> wireListHandle;
  std::vector<art::Ptr<recob::Wire>> wirelist;
  if (evt.getByLabel(fWireProducerLabel, wireListHandle)) {
    art::fill_ptr_vector(wirelist, wireListHandle);
  }

  if (rawdigitlist.empty() && wirelist.empty()) return;
  if (rawdigitlist.size() && wirelist.size()) return;

  // channel status
  lariov::ChannelStatusProvider const& channelStatus =
    art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  // ... Use the handle to get a particular (0th) element of collection.
  unsigned int dataSize;
  if (rawdigitlist.size()) {
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
    dataSize = digitVec0->Samples(); //size of raw data vectors
  }
  else {
    dataSize = (wirelist[0]->Signal()).size();
  }
  if (dataSize != detProp.ReadOutWindowSize()) {
    std::cout << "!!!!! Bad dataSize: " << dataSize << std::endl;
    return;
  }

  // ... Build a map from channel number -> rawdigitVec
  std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> rawdigitMap;
  raw::ChannelID_t chnum = raw::InvalidChannelID; // channel number
  if (rawdigitlist.size()) {
    for (size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter) {
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      chnum = digitVec->Channel();
      if (chnum == raw::InvalidChannelID) continue;
      rawdigitMap[chnum] = digitVec;
    }
  }
  // ... Build a map from channel number -> wire
  std::map<raw::ChannelID_t, art::Ptr<recob::Wire>> wireMap;
  if (wirelist.size()) {
    for (size_t ich = 0; ich < wirelist.size(); ++ich) {
      art::Ptr<recob::Wire> wire = wirelist[ich];
      chnum = wire->Channel();
      if (chnum == raw::InvalidChannelID) continue;
      wireMap[chnum] = wire;
    }
  }

  // ... Read in MC particle list
  art::Handle<std::vector<simb::MCParticle>> particleHandle;
  if (!evt.getByLabel(fSimulationProducerLabel, particleHandle)) {
    throw cet::exception("AnalysisExample")
      << " No simb::MCParticle objects in this event - "
      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  // ... Read in sim channel list
  auto simChannelHandle = evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);

  if (!simChannelHandle->size()) return;

  // ... Create a map of track IDs to generator labels
  //Get a list of generator names.
  //std::vector<art::Handle<std::vector<simb::MCTruth>>> mcHandles;
  //evt.getManyByType(mcHandles);
  auto mcHandles = evt.getMany<std::vector<simb::MCTruth>>();
  std::vector<std::pair<int, std::string>> track_id_to_label;

  for (auto const& mcHandle : mcHandles) {
    const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
    art::FindManyP<simb::MCParticle> findMCParts(mcHandle, evt, fSimulationProducerLabel);
    std::vector<art::Ptr<simb::MCParticle>> mcParts = findMCParts.at(0);
    for (const art::Ptr<simb::MCParticle> ptr : mcParts) {
      int track_id = ptr->TrackId();
      gf->add(track_id, sModuleLabel);
    }
  }

  std::string dummystr6 = "none  ";
  std::string dummystr7 = "none   ";

  if (fSaveSignal) {
    // .. create a channel number to trackid-wire signal info map
    std::map<raw::ChannelID_t, std::map<int, WireSigInfo>> Ch2TrkWSInfoMap;

    // .. create a track ID to vector of channel numbers (in w/c this track deposited energy) map
    std::map<int, std::vector<raw::ChannelID_t>> Trk2ChVecMap;

    // ... Loop over simChannels
    for (auto const& channel : (*simChannelHandle)) {

      // .. get simChannel channel number
      const raw::ChannelID_t ch1 = channel.Channel();
      if (ch1 == raw::InvalidChannelID) continue;
      if (geo::PlaneGeo::ViewName(fgeom->View(ch1)) != fPlaneToDump[0]) continue;

      bool selectThisChannel = false;

      // .. create a track ID to wire signal info map
      std::map<int, WireSigInfo> Trk2WSInfoMap;

      // ... Loop over all ticks with ionization energy deposited
      auto const& timeSlices = channel.TDCIDEMap();
      for (auto const& timeSlice : timeSlices) {

        auto const& energyDeposits = timeSlice.second;
        auto const tpctime = timeSlice.first;
        unsigned int tdctick = static_cast<unsigned int>(clockData.TPCTDC2Tick(double(tpctime)));
        if (tdctick < 0 || tdctick > (dataSize - 1)) continue;

        // ... Loop over all energy depositions in this tick
        for (auto const& energyDeposit : energyDeposits) {

          if (!energyDeposit.trackID) continue;
          int trkid = energyDeposit.trackID;
          simb::MCParticle particle = PIS->TrackIdToMotherParticle(trkid);
          //std::cout << energyDeposit.trackID << " " << trkid << " " << particle.TrackId() << std::endl;

          // .. ignore this energy deposition if incident particle energy below some threshold
          if (particle.E() < fMinParticleEnergyGeV) continue;

          int eve_id = PIS->TrackIdToEveTrackId(trkid);
          if (!eve_id) continue;
          std::string genlab = gf->get_gen(eve_id);

          if (Trk2WSInfoMap.find(trkid) == Trk2WSInfoMap.end()) {
            WireSigInfo wsinf;
            wsinf.pdgcode = particle.PdgCode();
            wsinf.genlab = genlab;
            wsinf.procid = particle.Process();
            wsinf.tdcmin = dataSize - 1;
            wsinf.tdcmax = 0;
            wsinf.edep = 0.;
            wsinf.numel = 0;
            Trk2WSInfoMap.insert(std::pair<int, WireSigInfo>(trkid, wsinf));
          }
          if (tdctick < Trk2WSInfoMap.at(trkid).tdcmin) Trk2WSInfoMap.at(trkid).tdcmin = tdctick;
          if (tdctick > Trk2WSInfoMap.at(trkid).tdcmax) Trk2WSInfoMap.at(trkid).tdcmax = tdctick;
          Trk2WSInfoMap.at(trkid).edep += energyDeposit.energy;
          Trk2WSInfoMap.at(trkid).numel += energyDeposit.numElectrons;
        }
      }

      if (!Trk2WSInfoMap.empty()) {
        for (std::pair<int, WireSigInfo> itmap : Trk2WSInfoMap) {
          if (fSelectGenLabel != "ANY") {
            if (itmap.second.genlab != fSelectGenLabel) continue;
          }
          if (fSelectProcID != "ANY") {
            if (itmap.second.procid != fSelectProcID) continue;
          }
          if (fSelectPDGCode != 0) {
            if (itmap.second.pdgcode != fSelectPDGCode) continue;
          }
          itmap.second.genlab.resize(6, ' ');
          itmap.second.procid.resize(7, ' ');
          if (itmap.second.numel >= fMinNumberOfElectrons &&
              itmap.second.edep >= fMinEnergyDepositedMeV) {
            if (fMaxNumberOfElectrons >= 0 && itmap.second.numel >= fMaxNumberOfElectrons) {
              continue;
            }
            else {
              int trkid = itmap.first;
              if (Trk2ChVecMap.find(trkid) == Trk2ChVecMap.end()) {
                std::vector<raw::ChannelID_t> chvec;
                Trk2ChVecMap.insert(std::pair<int, std::vector<raw::ChannelID_t>>(trkid, chvec));
              }
              Trk2ChVecMap.at(trkid).push_back(ch1);
              selectThisChannel = true;
            }
          }
        } // loop over Trk2WSinfoMap
        if (selectThisChannel) {
          Ch2TrkWSInfoMap.insert(
            std::pair<raw::ChannelID_t, std::map<int, WireSigInfo>>(ch1, Trk2WSInfoMap));
        }
      } // if Trk2WSInfoMap not empty

    } // loop over SimChannels

    std::set<raw::ChannelID_t> selected_channels;

    // ... Now write out the signal waveforms for each track
    if (!Trk2ChVecMap.empty()) {
      for (auto const& ittrk : Trk2ChVecMap) {
        int i = fRandFlat.fireInt(
          ittrk.second.size()); // randomly select one channel with a signal from this particle
        chnum = ittrk.second[i];

        if (not selected_channels.insert(chnum).second) { continue; }

        std::map<raw::ChannelID_t, std::map<int, WireSigInfo>>::iterator itchn;
        itchn = Ch2TrkWSInfoMap.find(chnum);
        if (itchn != Ch2TrkWSInfoMap.end()) {

          std::vector<short> adcvec(dataSize); // vector to hold zero-padded full waveform

          if (rawdigitlist.size()) {
            auto search = rawdigitMap.find(chnum);
            if (search == rawdigitMap.end()) continue;
            art::Ptr<raw::RawDigit> rawdig = (*search).second;
            std::vector<short> rawadc(dataSize); // vector to hold uncompressed adc values later
            raw::Uncompress(rawdig->ADCs(), rawadc, rawdig->GetPedestal(), rawdig->Compression());
            for (size_t j = 0; j < rawadc.size(); ++j) {
              adcvec[j] = rawadc[j] - rawdig->GetPedestal();
            }
          }
          else if (wirelist.size()) {
            auto search = wireMap.find(chnum);
            if (search == wireMap.end()) continue;
            art::Ptr<recob::Wire> wire = (*search).second;
            const auto& signal = wire->Signal();
            for (size_t j = 0; j < adcvec.size(); ++j) {
              adcvec[j] = signal[j];
            }
          }

          // .. write out info for each peak
          // a full waveform has at least one peak; the output will save up to 5 peaks (if there is
          // only 1 peak, will fill the other 4 with 0);
          // for fShortWaveformSize: only use the first peak's start_tick

          if (fUseFullWaveform) {

            c2numpy_uint32(&npywriter, evt.id().event());
            c2numpy_uint32(&npywriter, chnum);
            c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(chnum)).c_str());
            c2numpy_uint16(&npywriter, itchn->second.size()); // size of Trk2WSInfoMap, or #peaks
            unsigned int icnt = 0;
            for (auto const& it : itchn->second) {
              c2numpy_int32(&npywriter, it.first);                  // trackid
              c2numpy_int32(&npywriter, it.second.pdgcode);         // pdgcode
              c2numpy_string(&npywriter, it.second.genlab.c_str()); // genlab
              c2numpy_string(&npywriter, it.second.procid.c_str()); // procid
              c2numpy_float32(&npywriter, it.second.edep);          // edepo
              c2numpy_uint32(&npywriter, it.second.numel);          // numelec

              c2numpy_uint16(&npywriter, it.second.tdcmin); // stck1
              c2numpy_uint16(&npywriter, it.second.tdcmax); // stc2

              icnt++;
              if (icnt == 5) break;
            }

            // .. pad with 0's if number of peaks less than 5
            for (unsigned int i = icnt; i < 5; ++i) {
              c2numpy_int32(&npywriter, 0);
              c2numpy_int32(&npywriter, 0);
              c2numpy_string(&npywriter, dummystr6.c_str());
              c2numpy_string(&npywriter, dummystr7.c_str());
              c2numpy_float32(&npywriter, 0.);
              c2numpy_uint32(&npywriter, 0);
              c2numpy_uint16(&npywriter, 0);
              c2numpy_uint16(&npywriter, 0);
            }

            for (unsigned int itck = 0; itck < dataSize; ++itck) {
              c2numpy_int16(&npywriter, adcvec[itck]);
            }
          }
          else {

            // .. first loop to find largest signal
            double EDep = 0.;
            unsigned int TDCMin, TDCMax;
            bool foundmaxsig = false;
            for (auto& it : itchn->second) {
              if (it.second.edep > EDep && it.second.numel > 0) {
                EDep = it.second.edep;
                TDCMin = it.second.tdcmin;
                TDCMax = it.second.tdcmax;
                foundmaxsig = true;
              }
            }
            if (foundmaxsig) {
              int sigtdc1, sigtdc2, sighwid, sigfwid, sigtdcm;
              if (fPlaneToDump != fCollectionPlaneLabel) {
                sigtdc1 = TDCMin - 14 / 2;
                sigtdc2 = TDCMax + 3 * 14 / 2;
              }
              else {
                sigtdc1 = TDCMin - 32 / 2;
                sigtdc2 = TDCMax + 32 / 2;
              }
              sigfwid = sigtdc2 - sigtdc1;
              sighwid = sigfwid / 2;
              sigtdcm = sigtdc1 + sighwid;

              int start_tick = -1;
              int end_tick = -1;
              // .. set window edges to contain the largest signal
              if (sigfwid < (int)fShortWaveformSize) {
                // --> case 1: signal range fits within window
                int dt = fShortWaveformSize - sigfwid;
                start_tick = sigtdc1 - dt * fRandFlat.fire(0, 1);
              }
              else {
                // --> case 2: signal range larger than window
                int mrgn = fShortWaveformSize / 20;
                int dt = fShortWaveformSize - 2 * mrgn;
                start_tick = sigtdcm - mrgn - dt * fRandFlat.fire(0, 1);
              }
              if (start_tick < 0) start_tick = 0;
              end_tick = start_tick + fShortWaveformSize - 1;
              if (end_tick > int(dataSize - 1)) {
                end_tick = dataSize - 1;
                start_tick = end_tick - fShortWaveformSize + 1;
              }

              c2numpy_uint32(&npywriter, evt.id().event());
              c2numpy_uint32(&npywriter, chnum);
              c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(chnum)).c_str());

              // .. second loop to select only signals that are within the window

              int it_trk[5], it_pdg[5], it_nel[5];
              unsigned int stck_1[5], stck_2[5];
              std::string it_glb[5], it_prc[5];
              double it_edp[5];

              unsigned int icnt = 0;

              for (auto& it : itchn->second) {
                if ((it.second.tdcmin >= (unsigned int)start_tick &&
                     it.second.tdcmin < (unsigned int)end_tick) ||
                    (it.second.tdcmax > (unsigned int)start_tick &&
                     it.second.tdcmax <= (unsigned int)end_tick)) {

                  it_trk[icnt] = it.first;
                  it_pdg[icnt] = it.second.pdgcode;
                  it_glb[icnt] = it.second.genlab;
                  it_prc[icnt] = it.second.procid;
                  it_edp[icnt] = it.second.edep;
                  it_nel[icnt] = it.second.numel;

                  unsigned int mintdc = it.second.tdcmin;
                  unsigned int maxtdc = it.second.tdcmax;
                  if (mintdc < (unsigned int)start_tick) mintdc = start_tick;
                  if (maxtdc > (unsigned int)end_tick) maxtdc = end_tick;

                  stck_1[icnt] = mintdc - start_tick;
                  stck_2[icnt] = maxtdc - start_tick;

                  icnt++;
                  if (icnt == 5) break;
                }
              }

              c2numpy_uint16(&npywriter, icnt); // number of peaks

              for (unsigned int i = 0; i < icnt; ++i) {
                c2numpy_int32(&npywriter, it_trk[i]);          // trackid
                c2numpy_int32(&npywriter, it_pdg[i]);          // pdgcode
                c2numpy_string(&npywriter, it_glb[i].c_str()); // genlab
                c2numpy_string(&npywriter, it_prc[i].c_str()); // procid
                c2numpy_float32(&npywriter, it_edp[i]);        // edepo
                c2numpy_uint32(&npywriter, it_nel[i]);         // numelec
                c2numpy_uint16(&npywriter, stck_1[i]);         // stck1
                c2numpy_uint16(&npywriter, stck_2[i]);         // stck2
              }

              // .. pad with 0's if number of peaks less than 5
              for (unsigned int i = icnt; i < 5; ++i) {
                c2numpy_int32(&npywriter, 0);
                c2numpy_int32(&npywriter, 0);
                c2numpy_string(&npywriter, dummystr6.c_str());
                c2numpy_string(&npywriter, dummystr7.c_str());
                c2numpy_float32(&npywriter, 0.);
                c2numpy_uint32(&npywriter, 0);
                c2numpy_uint16(&npywriter, 0);
                c2numpy_uint16(&npywriter, 0);
              }

              for (unsigned int itck = start_tick; itck < (start_tick + fShortWaveformSize);
                   ++itck) {
                c2numpy_int16(&npywriter, adcvec[itck]);
              }

            } // foundmaxsig
          }
        }
      }
    }
  }
  else {
    //save noise
    int noisechancount = 0;
    std::map<raw::ChannelID_t, bool> signalMap;
    for (auto const& channel : (*simChannelHandle)) {
      signalMap[channel.Channel()] = true;
    }
    // .. create a vector for shuffling the wire channel indices
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    size_t nchan = (rawdigitlist.empty() ? wirelist.size() : rawdigitlist.size());
    std::vector<size_t> randigitmap;
    for (size_t i = 0; i < nchan; ++i)
      randigitmap.push_back(i);
    std::shuffle(randigitmap.begin(), randigitmap.end(), std::mt19937(seed));

    for (size_t rdIter = 0; rdIter < (rawdigitlist.empty() ? wirelist.size() : rawdigitlist.size());
         ++rdIter) {

      if (noisechancount == fMaxNoiseChannelsPerEvent) break;

      std::vector<short> adcvec(dataSize); // vector to wire adc values
      if (rawdigitlist.size()) {
        size_t ranIdx = randigitmap[rdIter];
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, ranIdx);
        if (signalMap[digitVec->Channel()]) continue;

        std::vector<short> rawadc(dataSize); // vector to hold uncompressed adc values later
        if (geo::PlaneGeo::ViewName(fgeom->View(digitVec->Channel())) != fPlaneToDump[0]) continue;
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
        for (size_t j = 0; j < rawadc.size(); ++j) {
          adcvec[j] = rawadc[j] - digitVec->GetPedestal();
        }
        c2numpy_uint32(&npywriter, evt.id().event());
        c2numpy_uint32(&npywriter, digitVec->Channel());
        c2numpy_string(&npywriter,
                       geo::PlaneGeo::ViewName(fgeom->View(digitVec->Channel())).c_str());
      }
      else if (wirelist.size()) {
        size_t ranIdx = randigitmap[rdIter];
        art::Ptr<recob::Wire> wire = wirelist[ranIdx];
        if (signalMap[wire->Channel()]) continue;
        if (channelStatus.IsBad(wire->Channel())) continue;
        if (geo::PlaneGeo::ViewName(fgeom->View(wire->Channel())) != fPlaneToDump[0]) continue;
        const auto& signal = wire->Signal();
        for (size_t j = 0; j < adcvec.size(); ++j) {
          adcvec[j] = signal[j];
        }
        c2numpy_uint32(&npywriter, evt.id().event());
        c2numpy_uint32(&npywriter, wire->Channel());
        c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(wire->Channel())).c_str());
      }

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
}
DEFINE_ART_MODULE(nnet::RawWaveformDump)
