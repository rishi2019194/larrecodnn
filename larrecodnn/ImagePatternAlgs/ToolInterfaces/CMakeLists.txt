cet_make_library(LIBRARY_NAME PointIdAlgorithm INTERFACE
  SOURCE IPointIdAlg.h
  LIBRARIES INTERFACE
  larreco::RecoAlg_ImagePatternAlgs_DataProvider
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  cetlib_except::cetlib_except
)

cet_write_plugin_builder(lar::PointIdAlgorithm art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larrecodnn::PointIdAlgorithm)

cet_make_library(LIBRARY_NAME WaveformRecognizer INTERFACE
  SOURCE IWaveformRecog.h
  LIBRARIES INTERFACE
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)

if (UNIX)
  target_link_libraries(WaveformRecognizer INTERFACE -lm)
endif()

cet_write_plugin_builder(lar::WaveformRecognizer art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larrecodnn::WaveformRecognizer)

cet_make_library(LIBRARY_NAME WireframeRecognizer INTERFACE
  SOURCE IWireframeRecog.h
  LIBRARIES INTERFACE
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)

if (UNIX)
  target_link_libraries(WireframeRecognizer INTERFACE -lm)
endif()

cet_write_plugin_builder(lar::WireframeRecognizer art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larrecodnn::WireframeRecognizer)



cet_make_library(LIBRARY_NAME WaveformDenoiser INTERFACE
  SOURCE IWaveformDenoise.h
  LIBRARIES INTERFACE
  canvas::canvas
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  cetlib_except::cetlib_except
)

if (UNIX)
  target_link_libraries(WaveformDenoiser INTERFACE -lm)
endif()

cet_write_plugin_builder(lar::WaveformDenoiser art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL larrecodnn::WaveformDenoiser)

install_headers()
