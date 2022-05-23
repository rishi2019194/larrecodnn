# ======================================================================
#  larrecodnn main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larrecodnn/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larrecodnn VERSION 09.11.04 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)
find_package(Threads)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package(Range-v3 REQUIRED PUBLIC)
find_package (Eigen3 3.3...3.3.9 REQUIRED NO_MODULE)
find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package(nusimdata REQUIRED PUBLIC)
find_package(nurandom REQUIRED PUBLIC)
find_package(nug4 REQUIRED PUBLIC)
find_package(nutools REQUIRED PUBLIC)
find_package(PostgreSQL REQUIRED PUBLIC)
find_package(ROOT COMPONENTS Core Minuit Minuit2 REQUIRED PUBLIC)

find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( canvas_root_io )
find_ups_root()
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( larreco )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

if(( DEFINED ENV{TENSORFLOW_DIR} ) OR ( DEFINED ENV{TRITON_DIR} ))
  find_package(Protobuf REQUIRED)
endif ()
# tensorflow is optional
if( DEFINED ENV{TENSORFLOW_DIR} )
  # find tensorflow library
  find_ups_product( tensorflow )
  cet_find_library(TENSORFLOW NAMES tensorflow PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
  cet_find_library(TENSORFLOW_CC NAMES tensorflow_cc PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
  cet_find_library(FOR_TENSORFLOW NAMES tensorflow_framework PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
endif ()
# trtis_clients is optional
if( DEFINED ENV{TRITON_DIR} )
  find_ups_product(triton)
  find_ups_product(grpc)
endif ()

if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  # find tensorflow library
  find_ups_product( protobuf )
  cet_find_library(PROTOBUF NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
endif()

include_directories( $ENV{EIGEN_INC} )

# source
add_subdirectory(larrecodnn)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
