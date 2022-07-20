#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/triton_utils.h"
#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/Span.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include <algorithm>
#include <experimental/iterator>
#include <iostream>
#include <sstream>
#include <string>

namespace triton_utils {

  template <typename C>
  std::string printColl(const C& coll, const std::string& delim) {
    if (coll.empty())
      return "";
    using std::begin;
    using std::end;
    std::stringstream msg;
    //avoid trailing delim
    std::copy(begin(coll), end(coll), std::experimental::make_ostream_joiner(msg, delim));
    return msg.str();
  }

  void throwIfError(const Error& err, std::string_view msg) {
    if (!err.IsOk())
      throw cet::exception("TritonServerFailure") << msg << ": " << err;
  }

  bool warnIfError(const Error& err, std::string_view msg) {
    if (!err.IsOk())
      MF_LOG_WARNING("TritonServerWarning") << msg << ": " << err;
    return err.IsOk();
  }

}  // namespace triton_utils

template std::string triton_utils::printColl(const triton_span::Span<std::vector<int64_t>::const_iterator>& coll,
                                             const std::string& delim);
template std::string triton_utils::printColl(const std::vector<uint8_t>& coll, const std::string& delim);
template std::string triton_utils::printColl(const std::vector<float>& coll, const std::string& delim);
template std::string triton_utils::printColl(const std::unordered_set<std::string>& coll, const std::string& delim);
