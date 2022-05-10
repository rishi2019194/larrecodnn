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

#include "dunereco/CVN/func/PixelMap.h"
#include "dunereco/CVN/func/SparsePixelMap.h"
#include "dunereco/CVN/func/Boundary.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace recob{class Hit, class Wire}; 
namespace sim{class SimChannel};

namespace cvn
{
  typedef std::vector<std::map<double, double>> Waveform;  
  
  class HitHelper
  {
  public:
    HitHelper(recob::Hit hit, double thresh = 0.): 
      fHit(hit), fThreshold(thresh) 
    {
      fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    }

    virtual Waveform GetWaveform();
    virtual geo::WireID GetID();
 
  private:
    recob::Hit fHit;
    double fThreshold;
    geo::GeometryCore const* fGeometry;
  };

  class WireHelper
  {
  public:
    WireHelper(recob::Wire wire, double thresh = 0.):
      fWire(wire), fThreshold(thresh)
    {
      fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    }

    virtual Waveform GetWaveform();
    virtual geo::WireID GetID();
  
  private:
    recob::Wire fWire;
    double fThreshold;
    geo::GeometryCore const* fGeometry;
  };

  class SimChannelHelper
  {
  public:
    SimChannelHelper(sim::SimChannel simchan, double thresh = 0.):
      fSimchan(simchan), fThreshold(thresh)
    {
      fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    }

    virtual Waveform GetWaveform();
    virtual geo::WireID GetID();
  
  private:
    sim::SimChannel fSimchan;
    double fThreshold;
    geo::GeometryCore const* fGeometry;
  };

  
  /// Producer algorithm for PixelMap, input to CVN neural net
  template <class T, class U> class PixelMapProducer
  {
  public:
    PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes, double threshold = 0.);
    PixelMapProducer();

    void SetMultipleDrifts() const {fMultipleDrifts = true;}
    
    /// Get boundaries for pixel map representation of cluster
    Boundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                            const std::vector< const T* >& cluster);

    virtual void ConvertLocaltoGlobal(geo::WireID wireid, 
                                      unsigned int &globalWire, unsigned int &globalPlane) const;    
    
    virtual void ConvertLocaltoGlobalTDC(geo::WireID wireid, double localTDC, 
                                         unsigned int &globalWire, unsigned int &globalPlane, 
                                         double &globalTDC) const;    

    unsigned int NWire() const {return fNWire;};
    unsigned int NTdc() const {return fNTdc;};
    double TRes() const {return fTRes;};

    PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                       const std::vector< art::Ptr< T > >& slice);
    PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                       const std::vector< const T* >& slice);

    PixelMap CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,
                                    const std::vector< const T* >& cluster,
                                    const Boundary& bound);

  private:
    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    double            fTRes;   ///< Timing resolution for pixel map
    unsigned int      fTotHits; ///< Total hits in the pixel map
    double            fThreshold; ///< Charge threshold to consider for hits/waveforms etc
    bool              fMultipleDrifts; ///< True if making the pixel map requires handling for multiple drift regions

    geo::GeometryCore const* fGeometry;
  };

  typedef PixelMapProducer<recob::Hit, HitHelper> PixelMapHitProducer;
  typedef PixelMapProducer<recob::Wire, WireHelper> PixelMapWireProducer;
  typedef PixelMapProducer<sim::SimChannel, SimChannelHelper> PixelMapSimProducer;

}

#endif  // CVN_PIXELMAPPRODUCER_H
