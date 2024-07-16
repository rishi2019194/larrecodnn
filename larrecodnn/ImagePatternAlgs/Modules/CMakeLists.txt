cet_make_library(LIBRARY_NAME EmTrack INTERFACE
  SOURCE
  EmTrack.h
  LIBRARIES INTERFACE
  larrecodnn::PointIdAlgorithm
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art_plugin_support::toolMaker
  art::Framework_Core
  art::Framework_Services_Registry
  art::Framework_Services_System_TriggerNamesService_service
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  cetlib::container_algorithms
  cetlib_except::cetlib_except
)

cet_build_plugin(CheckCNNScore art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::ArtDataHelper
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(EmTrackClusterId2outTl art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::EmTrack
  fhiclcpp::fhiclcpp
)

cet_build_plugin(EmTrackClusterId3outTl art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::EmTrack
  fhiclcpp::fhiclcpp
)

cet_build_plugin(EmTrackMichelIdTl art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::EmTrack
  fhiclcpp::fhiclcpp
)

cet_build_plugin(EvaluateROIEff art::EDAnalyzer
  LIBRARIES PRIVATE
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(PointIdEffTest art::EDAnalyzer
  LIBRARIES PRIVATE
  larreco::Calorimetry
  larsim::Simulation_LArG4Parameters_service
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  cetlib_except::cetlib_except
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(RawWaveformClnSigDump art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  CLHEP::Random
)

cet_build_plugin(RawWaveformDump art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  CLHEP::Random
)

cet_build_plugin(NoiseWaveformDump art::EDAnalyzer
  LIBRARIES PRIVATE
  larsim::MCCheater_ParticleInventoryService_service
  larevt::ChannelStatusProvider
  larevt::ChannelStatusService
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::Simulation
  larcoreobj::SimpleTypesAndConstants
  nusimdata::SimulationBase
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  CLHEP::Random
)

cet_build_plugin(WaveformRoiFinder art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::WaveformRecognizer
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(WaveformDenoiseTest art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::WaveformRecognizer
  larrecodnn::WaveformDenoiser
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardataobj::RawData
  lardataobj::RecoBase
  art_plugin_support::toolMaker
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

install_headers()
install_fhicl()
install_source()
