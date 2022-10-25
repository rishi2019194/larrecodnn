art_make(BASENAME_ONLY
#  LIBRARY_NAME      CVNArt
  LIB_LIBRARIES larrecodnn_CVN_func
  larrecodnn_CVN_tf
  larrecodnn_CVN_module_helpers
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  ART_PERSISTENCY_COMMON
  canvas::canvas
  art::Persistency_Provenance
  ART_UTILITIES
  messagefacility::MF_MessageLogger
  z
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  lardataobj_RecoBase
  lardata_Utilities
  nusimdata::SimulationBase
  larsim_MCCheater_ParticleInventoryService_service
  Boost::filesystem
  MODULE_LIBRARIES larrecodnn_CVN_func
  larrecodnn_CVN_tf
  larrecodnn_CVN_module_helpers
  # larrecodnn_CVN_modules
  stdc++fs
  )

install_headers()
install_fhicl()
install_source()
