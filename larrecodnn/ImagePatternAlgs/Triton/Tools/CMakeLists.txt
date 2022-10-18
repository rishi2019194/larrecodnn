include (lar::PointIdAlgorithm)

cet_build_plugin(PointIdAlgSonicTriton lar::PointIdAlgorithm
  LIBRARIES PRIVATE
  larrecodnn_ImagePatternAlgs_NuSonic_Triton
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
)

cet_build_plugin(PointIdAlgTriton lar::PointIdAlgorithm
  LIBRARIES PRIVATE
  TRITON::grpcclient
  gRPC::grpc # Should be a transitive (INTERFACE) dependency of TRITON::grpcclient
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

install_source()
