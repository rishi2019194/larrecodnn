
find_package(ZLIB REQUIRED)

add_subdirectory(func)
add_subdirectory(module_helpers)
add_subdirectory(modules)
if (TensorFlow_FOUND)
  add_subdirectory(tf)
endif ()
