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

install_headers()
install_source()
