
cet_build_plugin(NuGraphInference art::EDProducer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ROOT_MATH}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    lardataobj::RecoBase
    ${TORCH_LIBRARIES}
    TorchScatter::TorchScatter
    IMPL_TARGET_VAR NuGraphInference_module
)

find_path(DELAUNATOR_INC delaunator.hpp HINTS ENV DELAUNATOR_INC REQUIRED DOC "Header location for delaunator_cpp")
target_include_directories(${NuGraphInference_module} PRIVATE ${DELAUNATOR_INC})

cet_build_plugin(NuSliceHitsProducer art::EDProducer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    nusimdata::SimulationBase
    lardataobj::RecoBase
)

cet_build_plugin(NuGraphAnalyzer art::EDAnalyzer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    art_root_io::TFileService_service
    ROOT::Tree
    lardata::RecoBaseProxy
)



cet_build_plugin(NuGraphInferenceSonicTriton art::EDProducer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ROOT_MATH}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    lardataobj::RecoBase
    ${TORCH_LIBRARIES}
    TorchScatter::TorchScatter
    TritonClient::grpcclient
    gRPC::grpc # Should be a transitive (INTERFACE) dependency of TritonClient::grpcclient
    IMPL_TARGET_VAR NuGraphInferenceSonicTriton_module
)

cet_build_plugin(NuGraphInferenceTriton art::EDProducer
    LIBRARIES PRIVATE
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ROOT_MATH}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}
    lardataobj::RecoBase
    ${TORCH_LIBRARIES}
    TorchScatter::TorchScatter
    TritonClient::grpcclient
    gRPC::grpc # Should be a transitive (INTERFACE) dependency of TritonClient::grpcclient
    IMPL_TARGET_VAR NuGraphInferenceTriton_module
)

install_headers()
install_fhicl()
install_source()
