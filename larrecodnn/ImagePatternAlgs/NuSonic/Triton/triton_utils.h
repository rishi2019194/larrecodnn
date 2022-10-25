#ifndef NuSonic_Triton_triton_utils
#define NuSonic_Triton_triton_utils

#include <string>
#include <string_view>

#include "grpc_client.h"

namespace triton_utils {

  using Error = nvidia::inferenceserver::client::Error;

  template <typename C>
  std::string printColl(const C& coll, const std::string& delim = ", ");

  //helper to turn triton error into exception
  void throwIfError(const Error& err, std::string_view msg);

  //helper to turn triton error into warning
  bool warnIfError(const Error& err, std::string_view msg);

} // namespace triton_utils

#endif
