////////////////////////////////////////////////////////////////////////
/// \file    PixelMapProducer.h
/// \brief   PixelMapProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_PIXELMAPPRODUCER_H
#define CVN_PIXELMAPPRODUCER_H

#include <array>
#include <vector>

// Framework includes
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larrecodnn/CVN/func/Boundary.h"
#include "larrecodnn/CVN/func/PixelMap.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// namespace recob{class Hit, class Wire};
// namespace sim{class SimChannel};
//
namespace cvn {
  typedef std::vector<std::map<double, double>> Waveform;

  class HitHelper {
  public:
    HitHelper(recob::Hit hit, double thresh = 0.) : fHit(hit), fThreshold(thresh)
    {
      fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    }

    virtual Waveform GetWaveform();
    virtual geo::WireID GetID();

  protected:
    recob::Hit fHit;
    double fThreshold;
    geo::GeometryCore const* fGeometry;
  };

  class WireHelper {
  public:
    WireHelper(recob::Wire wire, double thresh = 0.) : fWire(wire), fThreshold(thresh)
    {
      fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    }

    virtual Waveform GetWaveform();
    virtual geo::WireID GetID();

  protected:
    recob::Wire fWire;
    double fThreshold;
    geo::GeometryCore const* fGeometry;
  };

  class SimChannelHelper {
  public:
    SimChannelHelper(sim::SimChannel simchan, double thresh = 0.)
      : fSimchan(simchan), fThreshold(thresh)
    {
      fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    }

    virtual Waveform GetWaveform();
    virtual geo::WireID GetID();

  protected:
    sim::SimChannel fSimchan;
    double fThreshold;
    geo::GeometryCore const* fGeometry;
  };

  /// Producer algorithm for PixelMap, input to CVN neural net
  template <class T, class U>
  class PixelMapProducer {
  public:
    PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes, double threshold = 0.);
    PixelMapProducer();

    // overload constructor for inputs from fcl
    PixelMapProducer(const fhicl::ParameterSet& pset);

    void SetMultipleDrifts() { fMultipleDrifts = true; }
    unsigned int NROI() { return fTotHits; };

    /// Get boundaries for pixel map representation of cluster
    virtual Boundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                    const std::vector<const T*>& cluster);

    virtual void ConvertLocaltoGlobal(geo::WireID wireid,
                                      unsigned int& globalWire,
                                      unsigned int& globalPlane) const;

    virtual void ConvertLocaltoGlobalTDC(geo::WireID wireid,
                                         double localTDC,
                                         unsigned int& globalWire,
                                         unsigned int& globalPlane,
                                         double& globalTDC) const;

    unsigned int NWire() const { return fNWire; }
    unsigned int NTdc() const { return fNTdc; }
    double TRes() const { return fTRes; }

    virtual PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                               const std::vector<art::Ptr<T>>& cluster);
    virtual PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                               const std::vector<const T*>& cluster);

    virtual PixelMap CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,
                                            const std::vector<const T*>& cluster,
                                            const Boundary& bound);

  protected:
    unsigned int fNWire;   ///< Number of wires, length for pixel maps
    unsigned int fNTdc;    ///< Number of tdcs, width of pixel map
    double fTRes;          ///< Timing resolution for pixel map
    unsigned int fTotHits; ///< Total hits in the pixel map
    double fThreshold;     ///< Charge threshold to consider for hits/waveforms etc
    bool
      fMultipleDrifts; ///< True if making the pixel map requires handling for multiple drift regions

    geo::GeometryCore const* fGeometry;
  };

  typedef PixelMapProducer<recob::Hit, cvn::HitHelper> PixelMapHitProducer;
  typedef PixelMapProducer<recob::Wire, cvn::WireHelper> PixelMapWireProducer;
  typedef PixelMapProducer<sim::SimChannel, cvn::SimChannelHelper> PixelMapSimProducer;

}

#endif // CVN_PIXELMAPPRODUCER_H
