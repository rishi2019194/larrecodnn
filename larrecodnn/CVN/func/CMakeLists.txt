
cet_make_library(SOURCE
  AssignLabels.cxx
  Boundary.cxx
  CVNImageUtils.cxx
  HitType.cxx
  InteractionType.cxx
  PixelMap.cxx
  Result.cxx
  TrainingData.cxx
  LIBRARIES
  PUBLIC
  nusimdata::SimulationBase
  art::Framework_Principal
  ROOT::Hist
  PRIVATE
  lardataobj::RecoBase
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  art::Framework_Services_Registry
  art::Persistency_Provenance
  art::Persistency_Common
  art::Utilities
  messagefacility::MF_MessageLogger
  canvas::canvas
)

art_dictionary(DICTIONARY_LIBRARIES larrecodnn::CVN_func)

install_headers()
install_source()

