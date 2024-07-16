////////////////////////////////////////////////////////////////////////
/// \file    PixelMapProducer.h
/// \brief   PixelMapProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
//
//  Modifications to allow unwrapped collection view
//   - Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <list>
#include <numeric>
#include <ostream>

#include "TH2D.h"
#include "TVector2.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace lcvn {

  Waveform HitHelper::GetWaveform()
  {

    Waveform ret;
    double pe = fHit.Integral();
    if (pe > fThreshold) ret.push_back(std::map<double, double>({{fHit.PeakTime(), pe}}));

    return ret;
  }

  geo::WireID HitHelper::GetID()
  {
    return fHit.WireID();
  }

  Waveform WireHelper::GetWaveform()
  {

    Waveform ret;
    auto ROIs = fWire.SignalROI();
    if (!(ROIs.get_ranges().size())) return ret;

    for (auto iROI = ROIs.begin_range(); iROI != ROIs.end_range(); ++iROI) {
      auto& ROI = *iROI;
      std::map<double, double> pulse;
      for (int tick = ROI.begin_index(); tick < (int)ROI.end_index(); tick++) {
        if (!(ROI[tick] > fThreshold)) continue;
        pulse.insert(std::pair<double, double>((double)tick, (double)ROI[tick]));
      }
      ret.push_back(pulse);
    }
    return ret;
  }

  geo::WireID WireHelper::GetID()
  {
    geo::WireID ret;
    std::vector<geo::WireID> wireids = fGeometry->ChannelToWire(fWire.Channel());
    if (!wireids.size()) return ret;
    ret = wireids[0];

    if (wireids.size() > 1) {
      for (auto iwire : wireids)
        if (iwire.Plane == fWire.View()) ret = iwire;
    }
    return ret;
  }

  Waveform SimChannelHelper::GetWaveform()
  {

    Waveform ret;
    auto& ROIs = fSimchan.TDCIDEMap();
    if (!(ROIs.size())) return ret;

    for (auto iROI = ROIs.begin(); iROI != ROIs.end(); ++iROI) {
      auto& ROI = *iROI;
      auto tick = ROI.first;
      double charge = 0.005 * fSimchan.Charge(tick);

      if (!(charge > fThreshold)) continue;
      ret.push_back(std::map<double, double>({{(double)tick, charge}}));
    }
    return ret;
  }

  geo::WireID SimChannelHelper::GetID()
  {
    geo::WireID ret;
    std::vector<geo::WireID> wireids = fGeometry->ChannelToWire(fSimchan.Channel());
    if (!wireids.size()) return ret;
    ret = wireids[0];
    return ret;
  }

  template <class T, class U>
  PixelMapProducer<T, U>::PixelMapProducer(unsigned int nWire,
                                           unsigned int nTdc,
                                           double tRes,
                                           double threshold)
    : fNWire(nWire), fNTdc(nTdc), fTRes(tRes), fThreshold(threshold), fMultipleDrifts(false)
  {

    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  }

  template <class T, class U>
  PixelMapProducer<T, U>::PixelMapProducer()
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  }

  template <class T, class U>
  PixelMapProducer<T, U>::PixelMapProducer(const fhicl::ParameterSet& pset)
    : fNWire(pset.get<unsigned int>("WireLength"))
    , fNTdc(pset.get<unsigned int>("TdcWidth"))
    , fTRes(pset.get<double>("TimeResolution"))
    , fThreshold(pset.get<double>("Threshold"))
    , fMultipleDrifts(pset.get<bool>("MultipleDrifts"))
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  }

  template <class T, class U>
  PixelMap PixelMapProducer<T, U>::CreateMap(detinfo::DetectorPropertiesData const& detProp,
                                             const std::vector<art::Ptr<T>>& cluster)
  {
    std::vector<const T*> newCluster;
    for (const art::Ptr<T> hit : cluster) {
      newCluster.push_back(hit.get());
    }
    return CreateMap(detProp, newCluster);
  }

  template <class T, class U>
  PixelMap PixelMapProducer<T, U>::CreateMap(detinfo::DetectorPropertiesData const& detProp,
                                             const std::vector<const T*>& cluster)
  {
    Boundary bound = DefineBoundary(detProp, cluster);
    return CreateMapGivenBoundary(detProp, cluster, bound);
  }

  template <class T, class U>
  PixelMap PixelMapProducer<T, U>::CreateMapGivenBoundary(
    detinfo::DetectorPropertiesData const& detProp,
    const std::vector<const T*>& cluster,
    const Boundary& bound)
  {

    PixelMap pm(fNWire, fNTdc, bound);

    for (size_t iHit = 0; iHit < cluster.size(); ++iHit) {

      U wraphit(*(cluster[iHit]), fThreshold);
      Waveform wf = wraphit.GetWaveform();
      geo::WireID wireid = wraphit.GetID();

      unsigned int tempWire = wireid.Wire;
      unsigned int tempPlane = wireid.Plane;

      if (!fMultipleDrifts) ConvertLocaltoGlobal(wireid, tempWire, tempPlane);

      for (auto& pulse : wf) {
        // Leigh: Simple modification to unwrap the collection view wire plane
        for (auto& i : pulse) {
          const double pe = i.second;
          double temptdc = i.first;
          if (fMultipleDrifts)
            ConvertLocaltoGlobalTDC(wireid, i.first, tempWire, tempPlane, temptdc);

          const unsigned int wire = tempWire;
          const unsigned int wirePlane = tempPlane;
          const double tdc = temptdc;

          pm.Add(wire, tdc, wirePlane, pe);
        }
      }
    }
    pm.SetTotHits(fTotHits);
    return pm;
  }

  template <class T, class U>
  std::ostream& operator<<(std::ostream& os, const PixelMapProducer<T, U>& p)
  {
    os << "PixelMapProducer: " << p.NTdc() << " tdcs X  " << p.NWire() << " wires";
    return os;
  }

  template <class T, class U>
  Boundary PixelMapProducer<T, U>::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                                  const std::vector<const T*>& cluster)
  {

    std::vector<double> tmin_0;
    std::vector<double> tmin_1;
    std::vector<double> tmin_2;

    std::vector<int> wire_0, bwire_0;
    std::vector<int> wire_1, bwire_1;
    std::vector<int> wire_2, bwire_2;

    std::vector<double> tsum = {0., 0., 0.};
    std::vector<double> tsize = {0., 0., 0.};

    for (size_t iHit = 0; iHit < cluster.size(); ++iHit) {
      U wraphit(*(cluster[iHit]), fThreshold);
      Waveform wf = wraphit.GetWaveform();
      geo::WireID wireid = wraphit.GetID();

      unsigned int tempWire = wireid.Wire;
      unsigned int tempPlane = wireid.Plane;

      if (!fMultipleDrifts) ConvertLocaltoGlobal(wireid, tempWire, tempPlane);

      for (auto& pulse : wf) {
        double min_tick = (double)INT_MAX;
        for (auto& i : pulse) {
          double temptdc = i.first;

          if (fMultipleDrifts)
            ConvertLocaltoGlobalTDC(wireid, i.first, tempWire, tempPlane, temptdc);

          if (temptdc < min_tick) min_tick = temptdc;

          tsum[tempPlane] += temptdc;
          tsize[tempPlane] += 1.;
        }

        if (!(pulse.empty())) {
          if (tempPlane == 0) {
            tmin_0.push_back(min_tick);
            wire_0.push_back(tempWire);
          }
          if (tempPlane == 1) {
            tmin_1.push_back(min_tick);
            wire_0.push_back(tempWire);
          }
          if (tempPlane == 2) {
            tmin_2.push_back(min_tick);
            wire_0.push_back(tempWire);
          }
        }
      } // end loop over pulses on single wire
    }   // end loop over struck wires

    double tmean_0 = tsum[0] / tsize[0];
    double tmean_1 = tsum[1] / tsize[1];
    double tmean_2 = tsum[2] / tsize[2];

    for (int i = 0; i < (int)wire_0.size(); i++) {
      if (std::abs(tmin_0[i] - tmean_0) < (double)fTRes) bwire_0.push_back(wire_0[i]);
    }
    for (int i = 0; i < (int)wire_1.size(); i++) {
      if (std::abs(tmin_1[i] - tmean_1) < (double)fTRes) bwire_1.push_back(wire_1[i]);
    }
    for (int i = 0; i < (int)wire_2.size(); i++) {
      if (std::abs(tmin_2[i] - tmean_2) < (double)fTRes) bwire_2.push_back(wire_2[i]);
    }

    std::cout << "Boundary wire vector sizes: " << bwire_0.size() << ", " << bwire_1.size() << ", "
              << bwire_2.size() << std::endl;

    int minwire_0 = 0;
    int minwire_1 = 0;
    int minwire_2 = 0;
    auto minwireelement_0 = std::min_element(bwire_0.begin(), bwire_0.end());
    auto minwireelement_1 = std::min_element(bwire_1.begin(), bwire_1.end());
    auto minwireelement_2 = std::min_element(bwire_2.begin(), bwire_2.end());

    if (bwire_0.size() > 0) {
      minwire_0 = *minwireelement_0 - 1;
      std::cout << "minwire 0: " << (*minwireelement_0 - 1) << std::endl;
    }
    if (bwire_1.size() > 0) {
      minwire_1 = *minwireelement_1 - 1;
      std::cout << "minwire 1: " << (*minwireelement_1 - 1) << std::endl;
    }
    if (bwire_2.size() > 0) {
      minwire_2 = *minwireelement_2 - 1;
      std::cout << "minwire 2: " << (*minwireelement_2 - 1) << std::endl;
    }

    fTotHits = bwire_0.size() + bwire_1.size() + bwire_2.size();

    Boundary bound(fNWire, fTRes, minwire_0, minwire_1, minwire_2, tmean_0, tmean_1, tmean_2);

    return bound;
  }

  template <class T, class U>
  void PixelMapProducer<T, U>::ConvertLocaltoGlobal(geo::WireID wireid,
                                                    unsigned int& globalWire,
                                                    unsigned int& globalPlane) const
  {
    globalWire = wireid.Wire;
    globalPlane = wireid.Plane;
  }

  template <class T, class U>
  void PixelMapProducer<T, U>::ConvertLocaltoGlobalTDC(geo::WireID wireid,
                                                       double localTDC,
                                                       unsigned int& globalWire,
                                                       unsigned int& globalPlane,
                                                       double& globalTDC) const
  {
    globalWire = wireid.Wire;
    globalPlane = wireid.Plane;
    globalTDC = localTDC;
  }

  template class PixelMapProducer<recob::Hit, lcvn::HitHelper>;
  template class PixelMapProducer<recob::Wire, lcvn::WireHelper>;
  template class PixelMapProducer<sim::SimChannel, lcvn::SimChannelHelper>;

} // namespace lcvn
