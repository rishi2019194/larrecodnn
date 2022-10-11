# ======================================================================
#  larrecodnn main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larrecodnn/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(larrecodnn VERSION 09.14.02 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

# needed because we call cet_make_plugin_builder
cet_cmake_module_directories(Modules BINARY)

# <<<<<<< HEAD
# # <<<<<<< HEAD
# find_ups_product( art )
# find_ups_product( art_root_io )
# find_ups_product( canvas_root_io )
# find_ups_root()
# find_ups_product( larcoreobj )
# find_ups_product( larcorealg )
# find_ups_product( larcore )
# find_ups_product( lardataobj )
# find_ups_product( lardataalg )
# find_ups_product( lardata )
# find_ups_product( larevt )
# find_ups_product( larsim )
# find_ups_product( larreco )
# # =======
# add_subdirectory(Modules)
# =======
find_package(larfinder REQUIRED EXPORT)
# >>>>>>> upstream/develop

find_package(larreco REQUIRED EXPORT)
find_package(larsim REQUIRED EXPORT)
find_package(larevt REQUIRED EXPORT)
find_package(lardata REQUIRED EXPORT)
find_package(larcore REQUIRED EXPORT)
find_package(lardataobj REQUIRED EXPORT)
find_package(larcorealg REQUIRED EXPORT)
find_package(larcoreobj REQUIRED EXPORT)

find_package(nusimdata REQUIRED EXPORT)
find_package(nurandom REQUIRED EXPORT)

find_package(art_root_io REQUIRED EXPORT)
find_package(art REQUIRED EXPORT)
find_package(canvas REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)

find_package(CLHEP COMPONENTS Random REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core Hist Tree REQUIRED EXPORT)
find_package(TBB REQUIRED EXPORT)
find_package(TRITON QUIET EXPORT)
find_package(TensorFlow 2.6.0 QUIET EXPORT)
find_package(Threads REQUIRED EXPORT)
# >>>>>>> upstream/develop

# macros for dictionary and simple_plugin
include(BuildPlugins)

# # <<<<<<< HEAD
# if(( DEFINED ENV{TENSORFLOW_DIR} ) OR ( DEFINED ENV{TRITON_DIR} ))
#   find_package(Protobuf REQUIRED)
# endif ()
# # tensorflow is optional
# if( DEFINED ENV{TENSORFLOW_DIR} )
#   # find tensorflow library
#   find_ups_product( tensorflow )
#   cet_find_library(TENSORFLOW NAMES tensorflow PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
#   cet_find_library(TENSORFLOW_CC NAMES tensorflow_cc PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
#   cet_find_library(FOR_TENSORFLOW NAMES tensorflow_framework PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
# endif ()
# # trtis_clients is optional
# if( DEFINED ENV{TRITON_DIR} )
#   find_ups_product(triton)
#   find_ups_product(grpc)
# endif ()
#
# if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
#   # find tensorflow library
#   find_ups_product( protobuf )
#   cet_find_library(PROTOBUF NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
# endif()
#
# include_directories( $ENV{EIGEN_INC} )
#
# # =======
# # >>>>>>> upstream/develop
# source
add_subdirectory(larrecodnn)

# ups - table and config files

# packaging utility
cet_cmake_config()
