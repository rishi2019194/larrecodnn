////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PointIdAlg
// Authors:     D.Stefan (Dorota.Stefan@ncbj.gov.pl),         from DUNE, CERN/NCBJ, since May 2016
//              R.Sulej (Robert.Sulej@cern.ch),               from DUNE, FNAL/NCBJ, since May 2016
//              P.Plonski,                                    from DUNE, WUT,       since May 2016
//
//
// Point Identification Algorithm
//
//      Run CNN or MLP trained to classify a point in 2D projection. Various features can be
//      recognized, depending on the net model/weights used.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PointIdAlg_h
#define PointIdAlg_h

// LArSoft includes
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

#include "larrecodnn/ImagePatternAlgs/Tensorflow/TF/tf_graph.h"
#include "larrecodnn/ImagePatternAlgs/Keras/keras_model.h"
#include "larreco/RecoAlg/ImagePatternAlgs/DataProvider/DataProviderAlg.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Framework includes
namespace art { class Event; }
namespace fhicl { class ParameterSet; }
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// ROOT & C++
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace nnet {
  class ModelInterface;
  class KerasModelInterface;
  class TfModelInterface;
  class PointIdAlg;
  class TrainingDataAlg;
}

/// Interface class for various classifier models. Now MLP (NetMaker) and CNN (Keras with
/// simple cpp interface) are supported. Will add interface to Protobuf as soon as Tensorflow
/// may be used from UPS.
class nnet::ModelInterface {
public:
  virtual ~ModelInterface() {}

  virtual std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) = 0;
  virtual std::vector<std::vector<float>> Run(
    std::vector<std::vector<std::vector<float>>> const& inps,
    int samples = -1);

protected:
  std::string findFile(const char* fileName) const;
};
// ------------------------------------------------------

class nnet::KerasModelInterface : public nnet::ModelInterface {
public:
  KerasModelInterface(const char* modelFileName);

  std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) override;

private:
  keras::KerasModel m; // network model
};
// ------------------------------------------------------

class nnet::TfModelInterface : public nnet::ModelInterface {
public:
  TfModelInterface(const char* modelFileName);

  std::vector<std::vector<float>> Run(std::vector<std::vector<std::vector<float>>> const& inps,
                                      int samples = -1) override;
  std::vector<float> Run(std::vector<std::vector<float>> const& inp2d) override;

private:
  std::unique_ptr<tf::Graph> g; // network graph
};
// ------------------------------------------------------

class nnet::PointIdAlg : public img::DataProviderAlg {
public:
  struct Config : public img::DataProviderAlg::Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<std::string> NNetModelFile{Name("NNetModelFile"),
                                           Comment("Neural net model to apply.")};
    fhicl::Sequence<std::string> NNetOutputs{Name("NNetOutputs"),
                                             Comment("Labels of the network outputs.")};
    fhicl::Atom<unsigned int> PatchSizeW{Name("PatchSizeW"), Comment("How many wires in patch.")};

    fhicl::Atom<unsigned int> PatchSizeD{Name("PatchSizeD"),
                                         Comment("How many downsampled ADC entries in patch")};
  };

  PointIdAlg(const fhicl::ParameterSet& pset) : PointIdAlg(fhicl::Table<Config>(pset, {})()) {}

  PointIdAlg(const Config& config);

  ~PointIdAlg() override;

  /// network output labels
  std::vector<std::string> const&
  outputLabels() const
  {
    return fNNetOutputs;
  }

  /// calculate single-value prediction (2-class probability) for [wire, drift] point
  float predictIdValue(unsigned int wire, float drift, size_t outIdx = 0) const;

  /// calculate multi-class probabilities for [wire, drift] point
  std::vector<float> predictIdVector(unsigned int wire, float drift) const;

  std::vector<std::vector<float>> predictIdVectors(
    std::vector<std::pair<unsigned int, float>> points) const;

  static std::vector<float> flattenData2D(std::vector<std::vector<float>> const& patch);

  std::vector<std::vector<float>> const&
  patchData2D() const
  {
    return fWireDriftPatch;
  }
  std::vector<float>
  patchData1D() const
  {
    return flattenData2D(fWireDriftPatch);
  } // flat vector made of the patch data, wire after wire

  bool isInsideFiducialRegion(unsigned int wire, float drift) const;

  /// test if wire/drift coordinates point to the current patch (so maybe the cnn output
  /// does not need to be recalculated)
  bool isCurrentPatch(unsigned int wire, float drift) const;

  /// test if two wire/drift coordinates point to the same patch
  bool isSamePatch(unsigned int wire1, float drift1, unsigned int wire2, float drift2) const;

private:
  std::string fNNetModelFilePath;
  std::vector<std::string> fNNetOutputs;
  nnet::ModelInterface* fNNet;

  mutable std::vector<std::vector<float>> fWireDriftPatch; // patch data around the identified point
  size_t fPatchSizeW, fPatchSizeD;

  mutable size_t fCurrentWireIdx, fCurrentScaledDrift;
  bool
  bufferPatch(size_t wire, float drift, std::vector<std::vector<float>>& patch) const
  {
    if (fDownscaleFullView) {
      size_t sd = (size_t)(drift / fDriftWindow);
      if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == sd))
        return true; // still within the current position

      fCurrentWireIdx = wire;
      fCurrentScaledDrift = sd;

      return patchFromDownsampledView(wire, drift, fPatchSizeW, fPatchSizeD, patch);
    }
    else {
      if ((fCurrentWireIdx == wire) && (fCurrentScaledDrift == drift))
        return true; // still within the current position

      fCurrentWireIdx = wire;
      fCurrentScaledDrift = drift;

      return patchFromOriginalView(wire, drift, fPatchSizeW, fPatchSizeD, patch);
    }
  }
  bool
  bufferPatch(size_t wire, float drift) const
  {
    return bufferPatch(wire, drift, fWireDriftPatch);
  }
  void resizePatch();

  void
  deleteNNet()
  {
    if (fNNet) delete fNNet;
    fNNet = 0;
  }
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class nnet::TrainingDataAlg : public img::DataProviderAlg {
public:
  enum EMask {
    kNone = 0,
    kPdgMask = 0x00000FFF,  // pdg code mask
    kTypeMask = 0x0000F000, // track type mask
    kVtxMask = 0xFFFF0000   // vertex flags
  };

  enum ETrkType {
    kDelta = 0x1000,  // delta electron
    kMichel = 0x2000, // Michel electron
    kPriEl = 0x4000,  // primary electron
    kPriMu = 0x8000   // primary muon
  };

  enum EVtxId {
    kNuNC = 0x0010000,
    kNuCC = 0x0020000,
    kNuPri = 0x0040000, // nu interaction type
    kNuE = 0x0100000,
    kNuMu = 0x0200000,
    kNuTau = 0x0400000,        // nu flavor
    kHadr = 0x1000000,         // hadronic inelastic scattering
    kPi0 = 0x2000000,          // pi0 produced in this vertex
    kDecay = 0x4000000,        // point of particle decay
    kConv = 0x8000000,         // gamma conversion
    kElectronEnd = 0x10000000, // clear end of an electron
    kElastic = 0x20000000,     // Elastic scattering
    kInelastic = 0x40000000    // Inelastic scattering
  };

  struct Config : public img::DataProviderAlg::Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> WireLabel{Name("WireLabel"), Comment("Tag of recob::Wire.")};

    fhicl::Atom<art::InputTag> HitLabel{Name("HitLabel"), Comment("Tag of recob::Hit.")};

    fhicl::Atom<art::InputTag> TrackLabel{Name("TrackLabel"), Comment("Tag of recob::Track.")};

    fhicl::Atom<art::InputTag> SimulationLabel{Name("SimulationLabel"),
                                               Comment("Tag of simulation producer.")};

    fhicl::Atom<art::InputTag> SimChannelLabel{Name("SimChannelLabel"),
                                               Comment("Tag of sim::SimChannel producer.")};

    fhicl::Atom<bool> SaveVtxFlags{Name("SaveVtxFlags"),
                                   Comment("Include (or not) vertex info in PDG map.")};

    fhicl::Atom<unsigned int> AdcDelayTicks{
      Name("AdcDelayTicks"),
      Comment("ADC pulse peak delay in ticks (non-zero for not deconvoluted waveforms).")};
  };

  TrainingDataAlg(const fhicl::ParameterSet& pset)
    : TrainingDataAlg(fhicl::Table<Config>(pset, {})())
  {}

  TrainingDataAlg(const Config& config);

  ~TrainingDataAlg() override;

  void reconfigure(const Config& config);

  bool
  saveSimInfo() const
  {
    return fSaveSimInfo;
  }

  bool setEventData(
    const art::Event& event, // collect & downscale ADC's, charge deposits, pdg labels
    detinfo::DetectorClocksData const& clockData,
    detinfo::DetectorPropertiesData const& detProp,
    unsigned int plane,
    unsigned int tpc,
    unsigned int cryo);

  bool setDataEventData(
    const art::Event& event, // collect & downscale ADC's, charge deposits, pdg labels
    detinfo::DetectorClocksData const& clockData,
    detinfo::DetectorPropertiesData const& detProp,
    unsigned int plane,
    unsigned int tpc,
    unsigned int cryo);

  bool findCrop(float max_e_cut,
                unsigned int& w0,
                unsigned int& w1,
                unsigned int& d0,
                unsigned int& d1) const;

  double
  getEdepTot() const
  {
    return fEdepTot;
  } // [GeV]
  std::vector<float> const&
  wireEdep(size_t widx) const
  {
    return fWireDriftEdep[widx];
  }
  std::vector<int> const&
  wirePdg(size_t widx) const
  {
    return fWireDriftPdg[widx];
  }

protected:
  img::DataProviderAlgView resizeView(detinfo::DetectorClocksData const& clock_data,
                                      detinfo::DetectorPropertiesData const& det_prop,
                                      size_t wires,
                                      size_t drifts) override;

private:
  struct WireDrift // used to find MCParticle start/end 2D projections
  {
    size_t Wire;
    int Drift;
    unsigned int TPC;
    unsigned int Cryo;
  };

  WireDrift getProjection(detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp,
                          const TLorentzVector& tvec,
                          unsigned int plane) const;

  bool setWireEdepsAndLabels(std::vector<float> const& edeps,
                             std::vector<int> const& pdgs,
                             size_t wireIdx);

  void collectVtxFlags(
    std::unordered_map<size_t, std::unordered_map<int, int>>& wireToDriftToVtxFlags,
    detinfo::DetectorClocksData const& clockData,
    detinfo::DetectorPropertiesData const& detProp,
    const std::unordered_map<int, const simb::MCParticle*>& particleMap,
    unsigned int plane) const;

  static float
  particleRange2(const simb::MCParticle& particle)
  {
    float dx = particle.EndX() - particle.Vx();
    float dy = particle.EndY() - particle.Vy();
    float dz = particle.EndZ() - particle.Vz();
    return dx * dx + dy * dy + dz * dz;
  }
  bool isElectronEnd(const simb::MCParticle& particle,
                     const std::unordered_map<int, const simb::MCParticle*>& particleMap) const;

  bool isMuonDecaying(const simb::MCParticle& particle,
                      const std::unordered_map<int, const simb::MCParticle*>& particleMap) const;

  double fEdepTot; // [GeV]
  std::vector<std::vector<float>> fWireDriftEdep;
  std::vector<std::vector<int>> fWireDriftPdg;

  art::InputTag fWireProducerLabel;
  art::InputTag fHitProducerLabel;
  art::InputTag fTrackModuleLabel;
  art::InputTag fSimulationProducerLabel;
  art::InputTag fSimChannelProducerLabel;
  bool fSaveVtxFlags;
  bool fSaveSimInfo;

  unsigned int fAdcDelay;

  std::vector<size_t> fEventsPerBin;
};
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

#endif
