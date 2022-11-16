

cet_build_plugin(TFNetHandler art::tool
  LIBRARIES PRIVATE
  larrecodnn_CVN_func
  larrecodnn_CVN_interfaces
  larrecodnn_ImagePatternAlgs_Tensorflow_TF
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  ROOT::Hist
)
