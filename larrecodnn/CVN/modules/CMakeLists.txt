
cet_build_plugin(LArNuCVNZlibMaker art::EDAnalyzer
  LIBRARIES PRIVATE
  larrecodnn::CVN_func
  larrecodnn::CVN_interfaces
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art_root_io::TFileService_service
  canvas::canvas
  Boost::filesystem
)

cet_build_plugin(LArCVNHitMapper art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::CVN_interfaces
  lardataobj::RecoBase
)

cet_build_plugin(LArCVNSimMapper art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::CVN_interfaces
  nusimdata::SimulationBase
)

cet_build_plugin(LArCVNWireMapper art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::CVN_interfaces
  lardataobj::RecoBase
)

install_headers()
install_fhicl()
install_source()
