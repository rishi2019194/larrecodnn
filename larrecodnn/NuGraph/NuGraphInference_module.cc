////////////////////////////////////////////////////////////////////////
// Class:       NuGraphInference
// Plugin Type: producer (Unknown Unknown)
// File:        NuGraphInference_module.cc
//
// Generated at Tue Nov 14 14:41:30 2023 by Giuseppe Cerati using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <array>
#include <limits>
#include <memory>

#include "delaunator.hpp"
#include <torch/script.h>

#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h" //this creates a conflict with torch script if included before it...

class NuGraphInference;

using anab::FeatureVector;
using anab::MVADescription;
using recob::Hit;
using recob::SpacePoint;
using std::array;
using std::vector;

namespace {
  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

  template <typename T, size_t N>
  void softmax(std::array<T, N>& arr)
  {
    T m = -std::numeric_limits<T>::max();
    for (size_t i = 0; i < arr.size(); i++) {
      if (arr[i] > m) { m = arr[i]; }
    }
    T sum = 0.0;
    for (size_t i = 0; i < arr.size(); i++) {
      sum += expf(arr[i] - m);
    }
    T offset = m + logf(sum);
    for (size_t i = 0; i < arr.size(); i++) {
      arr[i] = expf(arr[i] - offset);
    }
    return;
  }
}

class NuGraphInference : public art::EDProducer {
public:
  explicit NuGraphInference(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  NuGraphInference(NuGraphInference const&) = delete;
  NuGraphInference(NuGraphInference&&) = delete;
  NuGraphInference& operator=(NuGraphInference const&) = delete;
  NuGraphInference& operator=(NuGraphInference&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  vector<std::string> planes;
  art::InputTag hitInput;
  art::InputTag spsInput;
  size_t minHits;
  bool debug;
  vector<vector<float>> avgs;
  vector<vector<float>> devs;
  bool filterDecoder;
  bool semanticDecoder;
  bool vertexDecoder;
  torch::jit::script::Module model;
};

NuGraphInference::NuGraphInference(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , planes(p.get<vector<std::string>>("planes"))
  , hitInput(p.get<art::InputTag>("hitInput"))
  , spsInput(p.get<art::InputTag>("spsInput"))
  , minHits(p.get<size_t>("minHits"))
  , debug(p.get<bool>("debug"))
  , filterDecoder(p.get<bool>("filterDecoder"))
  , semanticDecoder(p.get<bool>("semanticDecoder"))
  , vertexDecoder(p.get<bool>("vertexDecoder"))
{

  for (size_t ip = 0; ip < planes.size(); ++ip) {
    avgs.push_back(p.get<vector<float>>("avgs_" + planes[ip]));
    devs.push_back(p.get<vector<float>>("devs_" + planes[ip]));
  }

  if (filterDecoder) { produces<vector<FeatureVector<1>>>("filter"); }
  //
  if (semanticDecoder) {
    produces<vector<FeatureVector<5>>>("semantic");
    produces<MVADescription<5>>("semantic");
  }
  //
  if (vertexDecoder) { produces<vector<recob::Vertex>>("vertex"); }

  cet::search_path sp("FW_SEARCH_PATH");
  model = torch::jit::load(sp.find_file(p.get<std::string>("modelFileName")));
}

void NuGraphInference::produce(art::Event& e)
{

  art::Handle<vector<Hit>> hitListHandle;
  vector<art::Ptr<Hit>> hitlist;
  if (e.getByLabel(hitInput, hitListHandle)) { art::fill_ptr_vector(hitlist, hitListHandle); }

  std::unique_ptr<vector<FeatureVector<1>>> filtcol(
    new vector<FeatureVector<1>>(hitlist.size(), FeatureVector<1>(std::array<float, 1>({-1.}))));

  std::unique_ptr<vector<FeatureVector<5>>> semtcol(new vector<FeatureVector<5>>(
    hitlist.size(), FeatureVector<5>(std::array<float, 5>({-1., -1., -1., -1., -1.}))));
  std::unique_ptr<MVADescription<5>> semtdes(
    new MVADescription<5>(hitListHandle.provenance()->moduleLabel(),
                          "semantic",
                          {"MIP", "HIP", "shower", "michel", "diffuse"}));

  std::unique_ptr<vector<recob::Vertex>> vertcol(new vector<recob::Vertex>());

  if (debug) std::cout << "Hits size=" << hitlist.size() << std::endl;
  if (hitlist.size() < minHits) {
    if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
    if (semanticDecoder) {
      e.put(std::move(semtcol), "semantic");
      e.put(std::move(semtdes), "semantic");
    }
    if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
    return;
  }

  vector<vector<float>> nodeft_bare(planes.size(), vector<float>());
  vector<vector<float>> nodeft(planes.size(), vector<float>());
  vector<vector<double>> coords(planes.size(), vector<double>());
  vector<vector<size_t>> idsmap(planes.size(), vector<size_t>());
  vector<size_t> idsmapRev(hitlist.size(), hitlist.size());
  for (auto h : hitlist) {
    idsmap[h->View()].push_back(h.key());
    idsmapRev[h.key()] = idsmap[h->View()].size() - 1;
    coords[h->View()].push_back(h->PeakTime() * 0.055);
    coords[h->View()].push_back(h->WireID().Wire * 0.3);
    nodeft[h->View()].push_back((h->WireID().Wire * 0.3 - avgs[h->View()][0]) / devs[h->View()][0]);
    nodeft[h->View()].push_back((h->PeakTime() * 0.055 - avgs[h->View()][1]) / devs[h->View()][1]);
    nodeft[h->View()].push_back((h->Integral() - avgs[h->View()][2]) / devs[h->View()][2]);
    nodeft[h->View()].push_back((h->RMS() - avgs[h->View()][3]) / devs[h->View()][3]);
    nodeft_bare[h->View()].push_back(h->WireID().Wire * 0.3);
    nodeft_bare[h->View()].push_back(h->PeakTime() * 0.055);
    nodeft_bare[h->View()].push_back(h->Integral());
    nodeft_bare[h->View()].push_back(h->RMS());
  }

  struct Edge {
    size_t n1;
    size_t n2;
    bool operator==(const Edge& other) const
    {
      if (this->n1 == other.n1 && this->n2 == other.n2)
        return true;
      else
        return false;
    };
  };
  vector<vector<Edge>> edge2d(planes.size(), vector<Edge>());
  for (size_t p = 0; p < planes.size(); p++) {
    if (debug) std::cout << "Plane " << p << " has N hits=" << coords[p].size() / 2 << std::endl;
    if (coords[p].size() / 2 < 3) continue;
    delaunator::Delaunator d(coords[p]);
    if (debug) std::cout << "Found N triangles=" << d.triangles.size() / 3 << std::endl;
    for (std::size_t i = 0; i < d.triangles.size(); i += 3) {
      //create edges in both directions
      Edge e;
      e.n1 = d.triangles[i];
      e.n2 = d.triangles[i + 1];
      edge2d[p].push_back(e);
      e.n1 = d.triangles[i + 1];
      e.n2 = d.triangles[i];
      edge2d[p].push_back(e);
      //
      e.n1 = d.triangles[i];
      e.n2 = d.triangles[i + 2];
      edge2d[p].push_back(e);
      e.n1 = d.triangles[i + 2];
      e.n2 = d.triangles[i];
      edge2d[p].push_back(e);
      //
      e.n1 = d.triangles[i + 1];
      e.n2 = d.triangles[i + 2];
      edge2d[p].push_back(e);
      e.n1 = d.triangles[i + 2];
      e.n2 = d.triangles[i + 1];
      edge2d[p].push_back(e);
      //
    }
    //sort and cleanup duplicate edges
    std::sort(edge2d[p].begin(), edge2d[p].end(), [](const auto& i, const auto& j) {
      return (i.n1 != j.n1 ? i.n1 < j.n1 : i.n2 < j.n2);
    });
    if (debug) {
      for (auto& e : edge2d[p]) {
        std::cout << "sorted plane=" << p << " e1=" << e.n1 << " e2=" << e.n2 << std::endl;
      }
    }
    edge2d[p].erase(std::unique(edge2d[p].begin(), edge2d[p].end()), edge2d[p].end());
  }

  if (debug) {
    for (size_t p = 0; p < planes.size(); p++) {
      for (auto& e : edge2d[p]) {
        std::cout << " plane=" << p << " e1=" << e.n1 << " e2=" << e.n2 << std::endl;
      }
    }
  }

  // Get spacepoints from the event record
  art::Handle<vector<SpacePoint>> spListHandle;
  vector<art::Ptr<SpacePoint>> splist;
  if (e.getByLabel(spsInput, spListHandle)) { art::fill_ptr_vector(splist, spListHandle); }
  // Get assocations from spacepoints to hits
  vector<vector<art::Ptr<Hit>>> sp2Hit(splist.size());
  if (splist.size() > 0) {
    art::FindManyP<Hit> fmp(spListHandle, e, "sps");
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    }
  }

  //Edges are the same as in pyg, but order is not identical.
  //It should not matter but better verify that output is indeed the same.
  vector<vector<Edge>> edge3d(planes.size(), vector<Edge>());
  for (size_t i = 0; i < splist.size(); ++i) {
    for (size_t j = 0; j < sp2Hit[i].size(); ++j) {
      Edge e;
      e.n1 = idsmapRev[sp2Hit[i][j].key()];
      e.n2 = i;
      edge3d[sp2Hit[i][j]->View()].push_back(e);
    }
  }

  auto x = torch::Dict<std::string, torch::Tensor>();
  auto batch = torch::Dict<std::string, torch::Tensor>();
  for (size_t p = 0; p < planes.size(); p++) {
    long int dim = nodeft[p].size() / 4;
    torch::Tensor ix = torch::zeros({dim, 4}, torch::dtype(torch::kFloat32));
    if (debug) {
      std::cout << "plane=" << p << std::endl;
      std::cout << std::fixed;
      std::cout << std::setprecision(4);
      std::cout << "before, plane=" << planes[p] << std::endl;
      for (size_t n = 0; n < nodeft_bare[p].size(); n = n + 4) {
        std::cout << nodeft_bare[p][n] << " " << nodeft_bare[p][n + 1] << " "
                  << nodeft_bare[p][n + 2] << " " << nodeft_bare[p][n + 3] << " " << std::endl;
      }
      std::cout << std::scientific;
      std::cout << "after, plane=" << planes[p] << std::endl;
      for (size_t n = 0; n < nodeft[p].size(); n = n + 4) {
        std::cout << nodeft[p][n] << " " << nodeft[p][n + 1] << " " << nodeft[p][n + 2] << " "
                  << nodeft[p][n + 3] << " " << std::endl;
      }
    }
    for (size_t n = 0; n < nodeft[p].size(); n = n + 4) {
      ix[n / 4][0] = nodeft[p][n];
      ix[n / 4][1] = nodeft[p][n + 1];
      ix[n / 4][2] = nodeft[p][n + 2];
      ix[n / 4][3] = nodeft[p][n + 3];
    }
    x.insert(planes[p], ix);
    torch::Tensor ib = torch::zeros({dim}, torch::dtype(torch::kInt64));
    batch.insert(planes[p], ib);
  }

  auto edge_index_plane = torch::Dict<std::string, torch::Tensor>();
  for (size_t p = 0; p < planes.size(); p++) {
    long int dim = edge2d[p].size();
    torch::Tensor ix = torch::zeros({2, dim}, torch::dtype(torch::kInt64));
    for (size_t n = 0; n < edge2d[p].size(); n++) {
      ix[0][n] = int(edge2d[p][n].n1);
      ix[1][n] = int(edge2d[p][n].n2);
    }
    edge_index_plane.insert(planes[p], ix);
    if (debug) {
      std::cout << "plane=" << p << std::endl;
      std::cout << "2d edge size=" << edge2d[p].size() << std::endl;
      for (size_t n = 0; n < edge2d[p].size(); n++) {
        std::cout << edge2d[p][n].n1 << " ";
      }
      std::cout << std::endl;
      for (size_t n = 0; n < edge2d[p].size(); n++) {
        std::cout << edge2d[p][n].n2 << " ";
      }
      std::cout << std::endl;
    }
  }

  auto edge_index_nexus = torch::Dict<std::string, torch::Tensor>();
  for (size_t p = 0; p < planes.size(); p++) {
    long int dim = edge3d[p].size();
    torch::Tensor ix = torch::zeros({2, dim}, torch::dtype(torch::kInt64));
    for (size_t n = 0; n < edge3d[p].size(); n++) {
      ix[0][n] = int(edge3d[p][n].n1);
      ix[1][n] = int(edge3d[p][n].n2);
    }
    edge_index_nexus.insert(planes[p], ix);
    if (debug) {
      std::cout << "plane=" << p << std::endl;
      std::cout << "3d edge size=" << edge3d[p].size() << std::endl;
      for (size_t n = 0; n < edge3d[p].size(); n++) {
        std::cout << edge3d[p][n].n1 << " ";
      }
      std::cout << std::endl;
      for (size_t n = 0; n < edge3d[p].size(); n++) {
        std::cout << edge3d[p][n].n2 << " ";
      }
      std::cout << std::endl;
    }
  }

  long int spdim = splist.size();
  auto nexus = torch::empty({spdim, 0}, torch::dtype(torch::kFloat32));

  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(x);
  inputs.push_back(edge_index_plane);
  inputs.push_back(edge_index_nexus);
  inputs.push_back(nexus);
  inputs.push_back(batch);
  if (debug) std::cout << "FORWARD!" << std::endl;
  auto outputs = model.forward(inputs).toGenericDict();
  if (debug) std::cout << "output =" << outputs << std::endl;
  if (semanticDecoder) {
    for (size_t p = 0; p < planes.size(); p++) {
      torch::Tensor s = outputs.at("x_semantic").toGenericDict().at(planes[p]).toTensor();
      for (int i = 0; i < s.sizes()[0]; ++i) {
        size_t idx = idsmap[p][i];
        std::array<float, 5> input({s[i][0].item<float>(),
                                    s[i][1].item<float>(),
                                    s[i][2].item<float>(),
                                    s[i][3].item<float>(),
                                    s[i][4].item<float>()});
        softmax(input);
        FeatureVector<5> semt = FeatureVector<5>(input);
        (*semtcol)[idx] = semt;
      }
      if (debug) {
        for (int j = 0; j < 5; j++) {
          std::cout << "x_semantic category=" << j << " : ";
          for (size_t p = 0; p < planes.size(); p++) {
            torch::Tensor s = outputs.at("x_semantic").toGenericDict().at(planes[p]).toTensor();
            for (int i = 0; i < s.sizes()[0]; ++i)
              std::cout << s[i][j].item<float>() << ", ";
          }
          std::cout << std::endl;
        }
      }
    }
  }
  if (filterDecoder) {
    for (size_t p = 0; p < planes.size(); p++) {
      torch::Tensor f = outputs.at("x_filter").toGenericDict().at(planes[p]).toTensor();
      for (int i = 0; i < f.numel(); ++i) {
        size_t idx = idsmap[p][i];
        std::array<float, 1> input({f[i].item<float>()});
        (*filtcol)[idx] = FeatureVector<1>(input);
      }
    }
    if (debug) {
      std::cout << "x_filter : ";
      for (size_t p = 0; p < planes.size(); p++) {
        torch::Tensor f = outputs.at("x_filter").toGenericDict().at(planes[p]).toTensor();
        for (int i = 0; i < f.numel(); ++i)
          std::cout << f[i].item<float>() << ", ";
      }
      std::cout << std::endl;
    }
  }
  if (vertexDecoder) {
    torch::Tensor v = outputs.at("x_vtx").toGenericDict().at("evt").toTensor()[0];
    double vpos[3];
    vpos[0] = v[0].item<float>();
    vpos[1] = v[1].item<float>();
    vpos[2] = v[2].item<float>();
    vertcol->push_back(recob::Vertex(vpos));
  }

  if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
  if (semanticDecoder) {
    e.put(std::move(semtcol), "semantic");
    e.put(std::move(semtdes), "semantic");
  }
  if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
}

DEFINE_ART_MODULE(NuGraphInference)
