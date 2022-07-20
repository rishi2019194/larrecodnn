# <<<<<<< HEAD
include_directories( $ENV{TENSORFLOW_INC}/absl )

art_make(
          LIB_LIBRARIES
			fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
			TENSORFLOW
			protobuf::libprotobuf
                        FOR_TENSORFLOW
      #rt #workaround for larsoft issue 17960
        )
# =======
# cet_make_library(SOURCE
#   tf_graph.cc
#   LIBRARIES PRIVATE
#   TensorFlow::cc
#   TensorFlow::framework
# )
# >>>>>>> upstream/develop

install_headers()
install_source()
