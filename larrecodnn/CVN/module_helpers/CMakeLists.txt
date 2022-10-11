#for classes_def.xml

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# if (DEFINED ENV{TENSORFLOW_DIR})
# # For some reason we apparently need to care if larreco is checked out.
# #if (larreco_not_in_ups)
# #set (LARRECO_RECOALG_IMAGEPATTERNALGS_TF larrecodnn_ImagePatternAlgs_Tensorflow_TF)
# #else (larreco_not_in_ups)
# #cet_find_library( LARRECO_RECOALG_IMAGEPATTERNALGS_TF NAMES larrecodnn_ImagePatternAlgs_Tensorflow_TF PATHS ENV LARRECO_LIB NO_DEFAULT_PATH)
# #endif (larreco_not_in_ups)
# else(DEFINED ENV{TENSORFLOW_DIR})
# set (EXCLUDE_TF TFNetHandler.cxx)
# endif (DEFINED ENV{TENSORFLOW_DIR})

art_make( BASENAME_ONLY
#  LIBRARY_NAME     CVNFunc
  # EXCLUDE ${EXCLUDE_TF}
  LIB_LIBRARIES larrecodnn_CVN_func
  larrecodnn_CVN_tf
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  cetlib::cetlib cetlib_except
  canvas::canvas
  Boost::filesystem            
  
  ROOT_BASIC_LIB_LIST
  DICT_LIBRARIES   lardataobj_RecoBase
  # MODULE_LIBRARIES larrecodnn_CVN_func
  larrecodnn_CVN_func
  larrecodnn_CVN_tf
  larrecodnn_CVN_module_helpers
  stdc++fs
  )

install_headers()
install_source()
