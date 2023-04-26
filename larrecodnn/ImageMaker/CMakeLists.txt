
# at this time, hep_hpc does not define target libraries
include_directories($ENV{HEP_HPC_INC})

cet_build_plugin(SaveImageH5 art::EDAnalyzer
  LIBRARIES PRIVATE
  hep_hpc_hdf5
  HDF5::HDF5
)

cet_build_plugin(SavePiMu art::tool
  LIBRARIES PRIVATE
  art::Framework_Core
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RecoBase
  nusimdata::SimulationBase
  hep_hpc_hdf5
  HDF5::HDF5
)

install_headers()
install_fhicl()
install_source()
