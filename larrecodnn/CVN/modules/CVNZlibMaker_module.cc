////////////////////////////////////////////////////////////////////////
// \file    CVNZlibMaker_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
//           - wrote the zlib code used in this module
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <cstdlib>

#include "boost/filesystem.hpp"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Data products
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
// #include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/TrainingData.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"
#include "larrecodnn/CVN/module_helpers/ICVNZlibMaker.h"

// Compression
#include "zlib.h"
#include "math.h"

#include "TH1.h"
#include "larcoreobj/SummaryData/POTSummary.h"

namespace cvn {

  class CVNZlibMaker : public ICVNZlibMaker {
  public:

    explicit CVNZlibMaker(fhicl::ParameterSet const& pset);
    ~CVNZlibMaker();

    void beginJob();
    void endSubRun(art::SubRun const &sr);
    void analyze(const art::Event& evt);
    void reconfigure(const fhicl::ParameterSet& pset);

    double SimpleOscProb(const simb::MCFlux& flux, const simb::MCNeutrino& nu) const;
  
  private:

    bool fIsVD;
    unsigned int fTopologyHitsCut;
    std::string fGenieGenModuleLabel;

    void write_files(TrainingNuData td, std::string evtid);

    TH1D* hPOT;
    double fPOT;
    int fRun;
    int fSubRun; 

  };

  //......................................................................
  CVNZlibMaker::CVNZlibMaker(fhicl::ParameterSet const& pset)
    : ICVNZlibMaker(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNZlibMaker::~CVNZlibMaker()
  {  }

  //......................................................................
  void CVNZlibMaker::reconfigure(const fhicl::ParameterSet& pset)
  {
    ICVNZlibMaker::reconfigure(pset);

    fIsVD = pset.get<bool>("IsVD");
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");
    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");

  }

  //......................................................................
  void CVNZlibMaker::endSubRun(const art::SubRun & sr){

    std::string fPOTModuleLabel = "generator";
    fRun = sr.run();
    fSubRun = sr.subRun();

    art::Handle< sumdata::POTSummary > potListHandle;
    if(sr.getByLabel(fPOTModuleLabel,potListHandle))
      fPOT = potListHandle->totpot;
    else
      fPOT = 0.;
    if(hPOT) hPOT->Fill(0.5, fPOT);
  }

  //......................................................................
  double CVNZlibMaker::SimpleOscProb(const simb::MCFlux& flux, const simb::MCNeutrino& nu) const {

    if( nu.CCNC() == 1 && nu.Nu().PdgCode() == flux.fntype) return 1;
    if( nu.CCNC() == 1 && nu.Nu().PdgCode() != flux.fntype) return 0;
    
    double E = nu.Nu().E();
    int flavAfter = nu.Nu().PdgCode();
    int flavBefore = flux.fntype;
    const double L      = 1284.9;    // km
    const double ldm    = 2.40e-3; // large delta m^2, m23
    const double ss2t13 = 0.1;     // sin^2(2theta_13)
    const double ss2t23 = 1.;      // maximal sin^2(2theta_23)
    const double ssth23 = 0.5;     // sin^2(theta_23) corresponding to above

    // Signal
    if(abs(flavAfter) == 12 && abs(flavBefore) == 14){
      // Nue appearance
      return ss2t13*ssth23*std::pow(sin(1.267*ldm*L/E), 2.);
    }
    if(abs(flavAfter) == 14 && abs(flavBefore) == 14){
      // CC mu disappearance
      return 1. - ss2t23*std::pow(sin(1.267*ldm*L/E), 2.);
    }

    // Background
    if(abs(flavAfter) == 12 && abs(flavBefore) == 12){
      // Beam nue
      return 1. - ss2t13*std::pow(sin(1.267*ldm*L/E), 2.);
    }
    if(abs(flavAfter) == 14 && abs(flavBefore) == 12){
      // CC mu appearance
      return ss2t13*ssth23*std::pow(sin(1.267*ldm*L/E), 2.);
    }
    if(abs(flavAfter) == 16 && abs(flavBefore) == 14){
      //numu to nutau CC appearance
      return (1.-ss2t13)*ss2t23*std::pow(sin(1.267*ldm*L/E), 2.);
    }
    if(abs(flavAfter) == 16 && abs(flavBefore) == 12){
      //nue to nutau CC appearance
      return ss2t13*(1.-ssth23)*std::pow(sin(1.267*ldm*L/E), 2.);
    }
    if(abs(flavAfter) == 16 && abs(flavBefore) == 16){
      //nutau to nutau CC disappearance
      return 1.-(ss2t23-ss2t23*ss2t13+ss2t13-ss2t13*ssth23)*std::pow(sin(1.267*ldm*L/E), 2.);
    }

    // Don't know what this is
    return 0;

  }
  
  //......................................................................
  void CVNZlibMaker::beginJob()
  {
    ICVNZlibMaker::beginJob();

    art::ServiceHandle<art::TFileService> tfs;
    hPOT = tfs->make<TH1D>("TotalPOT", "Total POT;; POT", 1, 0, 1);
  }

  //......................................................................
  void CVNZlibMaker::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    std::vector<art::Ptr<cvn::PixelMap>> pixelmaps;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto h_pixelmaps = evt.getHandle<std::vector<cvn::PixelMap>>(itag1);
    if (h_pixelmaps)
      art::fill_ptr_vector(pixelmaps, h_pixelmaps);

    // If no pixel maps, quit
    if (pixelmaps.size() == 0) return;

    InteractionType interaction = kOther;

    // MC information
    std::vector<art::Ptr<simb::MCTruth>> mctruth_list;
    auto h_mctruth = evt.getHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
    if (h_mctruth)
      art::fill_ptr_vector(mctruth_list, h_mctruth);

    art::Ptr<simb::MCTruth> mctruth = mctruth_list[0];
    simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();

    // Hard-coding event weight for now
    // Should probably fix this at some point
    double event_weight = 1.;

    art::PtrVector<simb::MCTruth> pv;
    pv.push_back(mctruth);
    art::FindManyP<simb::MCFlux> fmFlux(pv, evt, fGenieGenModuleLabel);

    if( fmFlux.isValid() ) {
      std::vector<art::Ptr<simb::MCFlux>> fluxes = fmFlux.at(0);
      if (!fluxes.empty()) {
        const simb::MCFlux& flux = *fluxes[0];
        event_weight = SimpleOscProb(flux, true_neutrino);
      }
    }

    AssignLabels labels;

    interaction = labels.GetInteractionType(true_neutrino);
    labels.GetTopology(mctruth, fTopologyHitsCut);

    // True lepton and neutrino energies
    float nu_energy = true_neutrino.Nu().E();
    float lep_energy = true_neutrino.Lepton().E();

    // Put a containment cut here
    bool fApplyFidVol = true;
    bool isFid = true;
    // If outside the fiducial volume don't waste any time filling other variables
    if(fApplyFidVol){
      // Get the interaction vertex from the end point of the neutrino. This is
      // because the start point of the lepton doesn't make sense for taus as they
      // are decayed by the generator and not GEANT
      TVector3 vtx = true_neutrino.Nu().EndPosition().Vect();
      if(fIsVD)
        isFid = (fabs(vtx.X())<300 && fabs(vtx.Y())<680 && vtx.Z()>40 && vtx.Z()<850); // vd
      else
        isFid = (fabs(vtx.X())<310 && fabs(vtx.Y())<550 && vtx.Z()>50 && vtx.Z()<1244); // hd
      if(!isFid) return;
    }

    TDNuInfo info;
    info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
    info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), 
                                labels.GetNPizeros(), labels.GetNNeutrons(), 
                                labels.GetTopologyType(), labels.GetTopologyTypeAlt());

    TrainingNuData train(interaction, *pixelmaps[0], info);
    

    std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_h"+std::to_string(time(0));
    this->write_files(train, evtid);
  }

  //......................................................................
  void CVNZlibMaker::write_files(TrainingNuData td, std::string evtid)
  {
    // cropped from 2880 x 500 to 500 x 500 here 
    std::vector<unsigned char> pixel_array(3 * fPlaneLimit * fTDCLimit);
   
    fImage.SetPixelMapSize(td.fPMap.NWire(), td.fPMap.NTdc());
    fImage.ConvertPixelMapToPixelArray(td.fPMap, pixel_array);

    ulong src_len = 3 * fPlaneLimit * fTDCLimit; // pixelArray length
    ulong dest_len = compressBound(src_len);     // calculate size of the compressed data               
    char* ostream = (char *) malloc(dest_len);  // allocate memory for the compressed data

    int res = compress((Bytef *) ostream, &dest_len, (Bytef *) &pixel_array[0], src_len);

    // Buffer error

    if (res == Z_BUF_ERROR)
      std::cout << "Buffer too small!" << std::endl;

    // Memory error
    else if (res ==  Z_MEM_ERROR)
      std::cout << "Not enough memory for compression!" << std::endl;

    // Compression ok 
    else {

      // Create output files 
      std::string image_file_name = out_dir + "/event_" + evtid + ".gz";
      std::string info_file_name = out_dir + "/event_" +  evtid + ".info";

      std::ofstream image_file (image_file_name, std::ofstream::binary);
      std::ofstream info_file  (info_file_name);

      if(image_file.is_open() && info_file.is_open()) {

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
    
    free(ostream);  // free allocated memory

  } // cvn::CVNZlibMaker::write_files

DEFINE_ART_MODULE(cvn::CVNZlibMaker)
} // namespace cvn
