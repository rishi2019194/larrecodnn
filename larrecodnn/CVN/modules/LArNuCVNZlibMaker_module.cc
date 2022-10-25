////////////////////////////////////////////////////////////////////////
// \file    LArNuCVNZlibMaker_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
//           - wrote the zlib code used in this module
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <cstdlib>
#include <iostream>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "boost/filesystem.hpp"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"

// Data products
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/TrainingData.h"
#include "larrecodnn/CVN/module_helpers/ICVNZlibMaker.h"

// Compression
#include "math.h"
#include "zlib.h"

#include "TH1.h"
#include "larcoreobj/SummaryData/POTSummary.h"

namespace cvn {

  class LArNuCVNZlibMaker : public cvn::ICVNZlibMaker {
  public:
    explicit LArNuCVNZlibMaker(fhicl::ParameterSet const& pset);
    ~LArNuCVNZlibMaker();

    void beginJob();
    void endSubRun(art::SubRun const& sr);
    void analyze(const art::Event& evt);
    void reconfigure(const fhicl::ParameterSet& pset);

  private:
    unsigned int fTopologyHitsCut;
    std::string fGenieGenModuleLabel;
    bool fApplyFidVol;
    std::vector<double> fFidMinCoords;
    std::vector<double> fFidMaxCoords;

    void write_files(TrainingNuData td, std::string evtid);

    TH1D* hPOT;
    double fPOT;
    int fRun;
    int fSubRun;
  };

  //......................................................................
  LArNuCVNZlibMaker::LArNuCVNZlibMaker(fhicl::ParameterSet const& pset) : ICVNZlibMaker(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  LArNuCVNZlibMaker::~LArNuCVNZlibMaker() {}

  //......................................................................
  void LArNuCVNZlibMaker::reconfigure(const fhicl::ParameterSet& pset)
  {
    ICVNZlibMaker::reconfigure(pset);

    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");
    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
    fApplyFidVol = pset.get<bool>("ApplyFidVol");
    fFidMinCoords = pset.get<std::vector<double>>("FidMinCoords");
    fFidMaxCoords = pset.get<std::vector<double>>("FidMaxCoords");
  }

  //......................................................................
  void LArNuCVNZlibMaker::endSubRun(const art::SubRun& sr)
  {

    std::string fPOTModuleLabel = "generator";
    fRun = sr.run();
    fSubRun = sr.subRun();

    art::Handle<sumdata::POTSummary> potListHandle;
    if (sr.getByLabel(fPOTModuleLabel, potListHandle))
      fPOT = potListHandle->totpot;
    else
      fPOT = 0.;
    if (hPOT) hPOT->Fill(0.5, fPOT);
  }

  //......................................................................
  void LArNuCVNZlibMaker::beginJob()
  {
    ICVNZlibMaker::beginJob();

    art::ServiceHandle<art::TFileService> tfs;
    hPOT = tfs->make<TH1D>("TotalPOT", "Total POT;; POT", 1, 0, 1);
  }

  //......................................................................
  void LArNuCVNZlibMaker::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    std::vector<art::Ptr<cvn::PixelMap>> pixelmaps;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto h_pixelmaps = evt.getHandle<std::vector<cvn::PixelMap>>(itag1);
    if (h_pixelmaps) art::fill_ptr_vector(pixelmaps, h_pixelmaps);

    // If no pixel maps, quit
    if (pixelmaps.size() == 0) return;

    InteractionType interaction = kOther;

    // MC information
    std::vector<art::Ptr<simb::MCTruth>> mctruth_list;
    auto h_mctruth = evt.getHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
    if (h_mctruth) art::fill_ptr_vector(mctruth_list, h_mctruth);

    art::Ptr<simb::MCTruth> mctruth = mctruth_list[0];
    simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();

    // Hard-coding event weight for now
    // Should probably fix this at some point
    double event_weight = 1.;

    AssignLabels labels;

    interaction = labels.GetInteractionType(true_neutrino);
    labels.GetTopology(mctruth, fTopologyHitsCut);

    // True lepton and neutrino energies
    float nu_energy = true_neutrino.Nu().E();
    float lep_energy = true_neutrino.Lepton().E();

    // Put a containment cut here
    // If outside the fiducial volume don't waste any time filling other variables
    if (fApplyFidVol) {
      // Get the interaction vertex from the end point of the neutrino. This is
      // because the start point of the lepton doesn't make sense for taus as they
      // are decayed by the generator and not GEANT
      TVector3 vtx = true_neutrino.Nu().EndPosition().Vect();
      bool isFid = (vtx.X() > fFidMinCoords[0] && vtx.X() < fFidMaxCoords[0]) &&
                   (vtx.Y() > fFidMinCoords[1] && vtx.Y() < fFidMaxCoords[1]) &&
                   (vtx.Z() > fFidMinCoords[2] && vtx.Z() < fFidMaxCoords[2]);
      if (!isFid) return;
    }

    TDNuInfo info;
    info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
    info.SetTopologyInformation(labels.GetPDG(),
                                labels.GetNProtons(),
                                labels.GetNPions(),
                                labels.GetNPizeros(),
                                labels.GetNNeutrons(),
                                labels.GetTopologyType(),
                                labels.GetTopologyTypeAlt());

    TrainingNuData train(interaction, *pixelmaps[0], info);

    std::string evtid = "r" + std::to_string(evt.run()) + "_s" + std::to_string(evt.subRun()) +
                        "_e" + std::to_string(evt.event()) + "_h" + std::to_string(time(0));
    this->write_files(train, evtid);
  }

  //......................................................................
  void LArNuCVNZlibMaker::write_files(TrainingNuData td, std::string evtid)
  {
    // cropped from 2880 x 500 to 500 x 500 here
    std::vector<unsigned char> pixel_array(3 * fPlaneLimit * fTDCLimit);

    fImage.SetPixelMapSize(td.fPMap.NWire(), td.fPMap.NTdc());
    fImage.ConvertPixelMapToPixelArray(td.fPMap, pixel_array);

    ulong src_len = 3 * fPlaneLimit * fTDCLimit; // pixelArray length
    ulong dest_len = compressBound(src_len);     // calculate size of the compressed data
    char* ostream = (char*)malloc(dest_len);     // allocate memory for the compressed data

    int res = compress((Bytef*)ostream, &dest_len, (Bytef*)&pixel_array[0], src_len);

    // Buffer error

    if (res == Z_BUF_ERROR) std::cout << "Buffer too small!" << std::endl;

    // Memory error
    else if (res == Z_MEM_ERROR)
      std::cout << "Not enough memory for compression!" << std::endl;

    // Compression ok
    else {

      // Create output files
      std::string image_file_name = out_dir + "/event_" + evtid + ".gz";
      std::string info_file_name = out_dir + "/event_" + evtid + ".info";

      std::ofstream image_file(image_file_name, std::ofstream::binary);
      std::ofstream info_file(info_file_name);

      if (image_file.is_open() && info_file.is_open()) {

        // Write compressed data to file

        image_file.write(ostream, dest_len);

        image_file.close(); // close file

        // Write records to file

        // Category

        info_file << td.fInt << std::endl;
        info_file << td.fInfo << std::endl;
        info_file << td.fPMap.GetTotHits() << std::endl;

        info_file.close(); // close file
      }
      else {

        if (image_file.is_open())
          image_file.close();
        else
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << image_file_name << "!" << std::endl;

        if (info_file.is_open())
          info_file.close();
        else
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << info_file_name << "!" << std::endl;
      }
    }

    free(ostream); // free allocated memory

  } // cvn::LArNuCVNZlibMaker::write_files

  DEFINE_ART_MODULE(cvn::LArNuCVNZlibMaker)
} // namespace cvn
