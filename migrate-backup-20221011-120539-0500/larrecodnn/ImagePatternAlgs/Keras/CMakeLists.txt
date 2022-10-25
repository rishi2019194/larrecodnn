cet_make_library(SOURCE
  keras_model.cc
  LIBRARIES PRIVATE
  TBB::tbb
)

add_subdirectory(Tools)

install_headers()
install_source()
