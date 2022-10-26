////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       Graph
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
////              P.Plonski,                      from DUNE, WUT, Sept. 2017
////              T.Cai (tejinc@yorku.ca)         from DUNE, YorkU, March 2022
////
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Graph_h
#define Graph_h

#include <memory>
#include <string>
#include <vector>

namespace tensorflow {
  class Session;
  class Tensor;
  struct SavedModelBundle;
}

namespace tf {

  class Graph {
  public:
    int n_inputs = 1;
    int n_outputs = 1;
    
    static std::unique_ptr<Graph> create(const char* graph_file_name,
                                         const std::vector<std::string>& outputs = {},
                                         bool use_bundle = false,
                                         int ninputs = 1,
                                         int noutputs = 1)
    {
      bool success;
      std::unique_ptr<Graph> ptr(new Graph(graph_file_name, outputs, success, use_bundle, ninputs, noutputs));
      if (success) { return ptr; }
      else {
        return nullptr;
      }
    }

    ~Graph();

    std::vector<float> run(const std::vector<std::vector<float>>& x);

    // process vector of 3D inputs, return vector of 1D outputs; use all inputs
    // if samples = -1, or only the specified number of first samples
    std::vector<std::vector<float>> run(
      const std::vector<std::vector<std::vector<std::vector<float>>>>& x,
      long long int samples = -1);
    std::vector<std::vector<float>> run(const tensorflow::Tensor& x);
   
    // use versions for multiple inputs -- needed for CVN 
    std::vector<std::vector<std::vector<float>>> runMulti(
      const std::vector<std::vector<std::vector<std::vector<float>>>>& x,
      long long int samples = -1);
    std::vector<std::vector<std::vector<float>>> runMulti(const std::vector<tensorflow::Tensor>& x);

  private:
    /// Not-throwing constructor.
    Graph(const char* graph_file_name,
          const std::vector<std::string>& outputs,
          bool& success,
          bool use_bundle = false,
          int ninputs = 1,
          int noutputs = 1);

    tensorflow::Session* fSession;
    bool fUseBundle;
    tensorflow::SavedModelBundle* fBundle;
    std::string fInputName;
    std::vector<std::string> fInputNames;
    std::vector<std::string> fOutputNames;
  };

} // namespace tf

#endif
