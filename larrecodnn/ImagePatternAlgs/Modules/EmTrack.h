#ifndef EMTRACK_H
#define EMTRACK_H

#include "lardata/ArtDataHelper/MVAWriter.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larrecodnn/ImagePatternAlgs/ToolInterfaces/IPointIdAlg.h"

#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib/container_algorithms.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace nnet {
  template <size_t N>
  class EmTrack {
  public:
    using key = std::tuple<unsigned int, unsigned int, unsigned int>;
    using cryo_tpc_view_keymap = std::map<key, std::vector<size_t>>;

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<PointIdAlgTools::IPointIdAlg::Config> PointIdAlg{Name("PointIdAlg")};
      fhicl::Atom<size_t> BatchSize{Name("BatchSize"),
                                    Comment("number of samples processed in one batch")};

      fhicl::Atom<art::InputTag> WireLabel{
        Name("WireLabel"),
        Comment("tag of deconvoluted ADC on wires (recob::Wire)")};

      fhicl::Atom<art::InputTag> HitModuleLabel{
        Name("HitModuleLabel"),
        Comment("tag of hits to be EM/track / Michel tagged")};

      fhicl::Atom<art::InputTag> ClusterModuleLabel{
        Name("ClusterModuleLabel"),
        Comment("tag of clusters to be used as a source of EM/track / Michel "
                "tagged new clusters (incl. single-hit clusters ) using "
                "accumulated results from hits")};

      fhicl::Atom<art::InputTag> TrackModuleLabel{
        Name("TrackModuleLabel"),
        Comment("tag of 3D tracks to be EM/track / Michel tagged using "
                "accumulated results from hits in the best 2D projection")};

      fhicl::Sequence<int> Views{Name("Views"),
                                 Comment("tag clusters in selected views only, "
                                         "or in all views if empty list")};
    };
    explicit EmTrack(Config const& c, std::string const& s, art::ProducesCollector& pc);
    void produce(art::Event& e);

  private:
    bool isViewSelected(int view) const;
    const size_t fBatchSize;
    std::unique_ptr<PointIdAlgTools::IPointIdAlg> fPointIdAlgTool;
    using writer = anab::MVAWriter<N>;
    writer fMVAWriter;
    const art::InputTag fWireProducerLabel;
    const art::InputTag fHitModuleLabel;
    const art::InputTag fClusterModuleLabel;
    const art::InputTag fTrackModuleLabel;
    const bool fDoClusters;
    const bool fDoTracks;
    const std::vector<int> fViews;
    const art::InputTag fNewClustersTag; // input tag for the clusters produced by this module
    void make_clusters(art::Event& evt,
                       std::vector<art::Ptr<recob::Hit>> const& hitPtrList,
                       std::vector<char> const& hitInFA,
                       EmTrack::cryo_tpc_view_keymap const& hitMap);
    void make_tracks(art::Event const& evt, std::vector<char> const& hitInFA);
    cryo_tpc_view_keymap create_hitmap(std::vector<art::Ptr<recob::Hit>> const& hitPtrList) const;
    std::vector<char> classify_hits(art::Event const& evt,
                                    EmTrack::cryo_tpc_view_keymap const& hitMap,
                                    std::vector<art::Ptr<recob::Hit>> const& hitPtrList);
  };

  template <size_t N>
  void EmTrack<N>::make_clusters(art::Event& evt,
                                 std::vector<art::Ptr<recob::Hit>> const& hitPtrList,
                                 std::vector<char> const& hitInFA,
                                 EmTrack::cryo_tpc_view_keymap const& hitMap)
  {
    // **************** prepare for new clusters ****************
    auto clusters = std::make_unique<std::vector<recob::Cluster>>();
    auto clu2hit = std::make_unique<art::Assns<recob::Cluster, recob::Hit>>();

    // ************** get and sort input clusters ***************
    auto cluListHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
    std::vector<art::Ptr<recob::Cluster>> cluPtrList;
    art::fill_ptr_vector(cluPtrList, cluListHandle);

    EmTrack::cryo_tpc_view_keymap cluMap;
    for (auto const& c : cluPtrList) {
      unsigned int view = c->Plane().Plane;
      if (!isViewSelected(view)) continue;

      unsigned int cryo = c->Plane().Cryostat;
      unsigned int tpc = c->Plane().TPC;

      cluMap[{cryo, tpc, view}].push_back(c.key());
    }

    auto cluID = fMVAWriter.template initOutputs<recob::Cluster>(fNewClustersTag,
                                                                 fPointIdAlgTool->outputLabels());

    unsigned int cidx = 0; // new clusters index
    art::FindManyP<recob::Hit> hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
    std::vector<bool> hitUsed(hitPtrList.size(),
                              false);      // tag hits used in clusters
                                           // clang-format off
        for (auto const & [key, clusters_keys] : cluMap)
        {
            auto const& [cryo, tpc, view]= key;
                                           // clang-format on
      if (!isViewSelected(view)) continue; // should not happen, clusters were pre-selected

      for (size_t c : clusters_keys) // c is the Ptr< recob::Cluster >::key()
      {
        auto v = hitsFromClusters.at(c);
        if (v.empty()) continue;

        for (auto const& hit : v) {
          if (hitUsed[hit.key()]) {
            mf::LogWarning("EmTrack") << "hit already used in another cluster";
          }
          hitUsed[hit.key()] = true;
        }

        auto vout = fMVAWriter.template getOutput<recob::Hit>(
          v, [&](art::Ptr<recob::Hit> const& ptr) { return (float)hitInFA[ptr.key()]; });

        float pvalue = vout[0] / (vout[0] + vout[1]);
        mf::LogVerbatim("EmTrack") << "cluster in tpc:" << tpc << " view:" << view
                                   << " size:" << v.size() << " p:" << pvalue;

        clusters->emplace_back(recob::Cluster(0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              v.size(),
                                              0.0F,
                                              0.0F,
                                              cidx,
                                              (geo::View_t)view,
                                              v.front()->WireID().planeID()));
        util::CreateAssn(evt, *clusters, v, *clu2hit);
        cidx++;

        fMVAWriter.addOutput(cluID,
                             vout); // add copy of the input cluster
      }

      // (2b) make single-hit clusters
      // --------------------------------------------
      for (size_t h : hitMap.at({cryo, tpc, view})) // h is the Ptr< recob::Hit >::key()
      {
        if (hitUsed[h]) continue;

        auto vout = fMVAWriter.template getOutput<recob::Hit>(h);
        float pvalue = vout[0] / (vout[0] + vout[1]);

        mf::LogVerbatim("EmTrack") << "single hit in tpc:" << tpc << " view:" << view
                                   << " wire:" << hitPtrList[h]->WireID().Wire
                                   << " drift:" << hitPtrList[h]->PeakTime() << " p:" << pvalue;

        art::PtrVector<recob::Hit> cluster_hits;
        cluster_hits.push_back(hitPtrList[h]);
        clusters->emplace_back(recob::Cluster(0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              0.0F,
                                              1,
                                              0.0F,
                                              0.0F,
                                              cidx,
                                              (geo::View_t)view,
                                              hitPtrList[h]->WireID().planeID()));
        util::CreateAssn(evt, *clusters, cluster_hits, *clu2hit);
        cidx++;

        fMVAWriter.addOutput(cluID, vout); // add single-hit cluster tagging unclutered hit
      }
      mf::LogVerbatim("EmTrack") << "...produced " << cidx - clusters_keys.size()
                                 << " single-hit clusters.";
    }

    evt.put(std::move(clusters));
    evt.put(std::move(clu2hit));
  }

  /// make tracks
  template <size_t N>
  void EmTrack<N>::make_tracks(art::Event const& evt, std::vector<char> const& hitInFA)
  {
    auto trkListHandle = evt.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    art::FindManyP<recob::Hit> hitsFromTracks(trkListHandle, evt, fTrackModuleLabel);
    std::vector<std::vector<art::Ptr<recob::Hit>>> trkHitPtrList(trkListHandle->size());
    for (size_t t = 0; t < trkListHandle->size(); ++t) {
      auto v = hitsFromTracks.at(t);
      size_t nh[3] = {0, 0, 0};
      for (auto const& hptr : v) {
        ++nh[hptr->View()];
      }
      size_t best_view = 2;                                       // collection
      if ((nh[0] >= nh[1]) && (nh[0] > 2 * nh[2])) best_view = 0; // ind1
      if ((nh[1] >= nh[0]) && (nh[1] > 2 * nh[2])) best_view = 1; // ind2

      size_t k = 0;
      while (!isViewSelected(best_view)) {
        best_view = (best_view + 1) % 3;
        if (++k > 3) {
          throw cet::exception("EmTrack") << "No views selected at all?" << std::endl;
        }
      }

      for (auto const& hptr : v) {
        if (hptr->View() == best_view) trkHitPtrList[t].emplace_back(hptr);
      }
    }

    auto trkID = fMVAWriter.template initOutputs<recob::Track>(
      fTrackModuleLabel, trkHitPtrList.size(), fPointIdAlgTool->outputLabels());
    for (size_t t = 0; t < trkHitPtrList.size(); ++t) // t is the Ptr< recob::Track >::key()
    {
      auto vout = fMVAWriter.template getOutput<recob::Hit>(
        trkHitPtrList[t],
        [&](art::Ptr<recob::Hit> const& ptr) { return (float)hitInFA[ptr.key()]; });
      fMVAWriter.setOutput(trkID, t, vout);
    }
  }
  template <size_t N>
  typename EmTrack<N>::cryo_tpc_view_keymap EmTrack<N>::create_hitmap(
    std::vector<art::Ptr<recob::Hit>> const& hitPtrList) const
  {
    cryo_tpc_view_keymap hitMap;
    for (auto const& hptr : hitPtrList) {
      auto const& h = *hptr;
      unsigned int view = h.WireID().Plane;
      if (!isViewSelected(view)) continue;

      unsigned int cryo = h.WireID().Cryostat;
      unsigned int tpc = h.WireID().TPC;

      hitMap[{cryo, tpc, view}].push_back(hptr.key());
    }
    return hitMap;
  }

  template <size_t N>
  std::vector<char> EmTrack<N>::classify_hits(art::Event const& evt,
                                              EmTrack::cryo_tpc_view_keymap const& hitMap,
                                              std::vector<art::Ptr<recob::Hit>> const& hitPtrList)
  {
    auto hitID = fMVAWriter.template initOutputs<recob::Hit>(
      fHitModuleLabel, hitPtrList.size(), fPointIdAlgTool->outputLabels());

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    auto wireHandle = evt.getValidHandle<std::vector<recob::Wire>>(fWireProducerLabel);
    std::vector<char> hitInFA(hitPtrList.size(),
                              0); // tag hits in fid. area as 1, use 0 for hits
                                  // close to the projectrion edges
    for (auto const& [key, hits] : hitMap) {
      auto const& [cryo, tpc, view] = key;
      if (!isViewSelected(view)) continue; // should not happen, hits were selected

      fPointIdAlgTool->setWireDriftData(clockData, detProp, *wireHandle, view, tpc, cryo);

      // (1) do all hits in this plane
      // ------------------------------------------------
      for (size_t idx = 0; idx < hits.size(); idx += fBatchSize) {
        std::vector<std::pair<unsigned int, float>> points;
        std::vector<size_t> keys;
        for (size_t k = 0; k < fBatchSize; ++k) {
          if (idx + k >= hits.size()) { break; } // careful about the tail

          size_t h = hits[idx + k]; // h is the Ptr< recob::Hit >::key()
          const recob::Hit& hit = *(hitPtrList[h]);
          points.emplace_back(hit.WireID().Wire, hit.PeakTime());
          keys.push_back(h);
        }

        auto batch_out = fPointIdAlgTool->predictIdVectors(points);
        if (points.size() != batch_out.size()) {
          throw cet::exception("EmTrack") << "hits processing failed" << std::endl;
        }

        for (size_t k = 0; k < points.size(); ++k) {
          size_t h = keys[k];
          fMVAWriter.setOutput(hitID, h, batch_out[k]);
          if (fPointIdAlgTool->isInsideFiducialRegion(points[k].first, points[k].second)) {
            hitInFA[h] = 1;
          }
        }
      } // hits done
        // ------------------------------------------------------------------
    }
    return hitInFA;
  }
  // make sure fMVAWriter is getting a variable string
  template <size_t N>
  EmTrack<N>::EmTrack(EmTrack::Config const& config,
                      std::string const& module_label,
                      art::ProducesCollector& collector)
    : fBatchSize(config.BatchSize())
    , fPointIdAlgTool(art::make_tool<PointIdAlgTools::IPointIdAlg>(config.PointIdAlg.get_PSet()))
    , fMVAWriter(collector, "emtrkmichel")
    , fWireProducerLabel(config.WireLabel())
    , fHitModuleLabel(config.HitModuleLabel())
    , fClusterModuleLabel(config.ClusterModuleLabel())
    , fTrackModuleLabel(config.TrackModuleLabel())
    , fDoClusters(!fClusterModuleLabel.label().empty())
    , fDoTracks(!fTrackModuleLabel.label().empty())
    , fViews(config.Views())
    , fNewClustersTag(module_label,
                      "",
                      art::ServiceHandle<art::TriggerNamesService const>()->getProcessName())
  {
    fMVAWriter.template produces_using<recob::Hit>();

    if (!fClusterModuleLabel.label().empty()) {
      collector.produces<std::vector<recob::Cluster>>();
      collector.produces<art::Assns<recob::Cluster, recob::Hit>>();

      fMVAWriter.template produces_using<recob::Cluster>();
    }

    if (!fTrackModuleLabel.label().empty()) { fMVAWriter.template produces_using<recob::Track>(); }
  }
  // ------------------------------------------------------

  template <size_t N>
  void EmTrack<N>::produce(art::Event& evt)
  {
    mf::LogVerbatim("EmTrack") << "next event: " << evt.run() << " / " << evt.id().event();
    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);
    std::vector<art::Ptr<recob::Hit>> hitPtrList;
    art::fill_ptr_vector(hitPtrList, hitListHandle);
    const EmTrack::cryo_tpc_view_keymap hitMap = create_hitmap(hitPtrList);
    const std::vector<char> hitInFA = classify_hits(evt, hitMap, hitPtrList);

    if (fDoClusters) make_clusters(evt, hitPtrList, hitInFA, hitMap);

    if (fDoTracks) make_tracks(evt, hitInFA);
    fMVAWriter.saveOutputs(evt);
  }
  // ------------------------------------------------------

  template <size_t N>
  bool EmTrack<N>::isViewSelected(int view) const
  {
    if (fViews.empty()) return true;
    return cet::search_all(fViews, view);
  }
  // ------------------------------------------------------

}
#endif
