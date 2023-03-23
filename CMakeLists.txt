# ======================================================================
#  larrecodnn main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larrecodnn/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(larrecodnn VERSION 09.17.13 LANGUAGES CXX)

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

find_package(larfinder REQUIRED EXPORT)

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

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BuildPlugins)

# source
add_subdirectory(larrecodnn)

# ups - table and config files

# packaging utility
cet_cmake_config()
