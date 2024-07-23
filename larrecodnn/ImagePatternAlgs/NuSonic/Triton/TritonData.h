#ifndef NuSonic_Triton_TritonData
#define NuSonic_Triton_TritonData

#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/Span.h"
#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/triton_utils.h"

#include <algorithm>
#include <any>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include "grpc_client.h"
#include "triton/common/model_config.h"

namespace nic = triton::client;

namespace lartriton {

  //forward declaration
  class TritonClient;

  //aliases for local input and output types
  template <typename DT>
  using TritonInput = std::vector<std::vector<DT>>;
  template <typename DT>
  using TritonOutput = std::vector<triton_span::Span<const DT*>>;

  //store all the info needed for triton input and output
  template <typename IO>
  class TritonData {
  public:
    using Result = nic::InferResult;
    using TensorMetadata = inference::ModelMetadataResponse_TensorMetadata;
    using ShapeType = std::vector<int64_t>;
    using ShapeView = triton_span::Span<ShapeType::const_iterator>;

    //constructor
    TritonData(const std::string& name, const TensorMetadata& model_info, bool noBatch);

    //some members can be modified
    bool setShape(const ShapeType& newShape) { return setShape(newShape, true); }
    bool setShape(unsigned loc, int64_t val) { return setShape(loc, val, true); }

    //io accessors
    template <typename DT>
    void toServer(std::shared_ptr<TritonInput<DT>> ptr) {
        // Access the actual data inside the shared_ptr
        const auto& data_in = *ptr;

        // Assume 'data_' is a member variable or accessible globally
        // Set the shape of the data if necessary
        data_->SetShape({static_cast<long int>(data_in[0].size())});

        // Check byte size consistency
        if (byteSize_ != sizeof(DT)) {
            throw std::runtime_error("Byte size inconsistency");
        }

        // Iterate through each element in the batch
        for (unsigned i0 = 0; i0 < batchSize_; ++i0) {
            const DT* arr = data_in[i0].data(); // Access raw data pointer
            // Append raw data to Triton server
            triton_utils::throwIfError(
                data_->AppendRaw(reinterpret_cast<const uint8_t*>(arr), data_in[i0].size() * byteSize_),
                "Failed to set data for batch entry " + std::to_string(i0));
        }

        // Store the shared_ptr to keep it in scope
        holder_ = std::move(ptr);
    }

    template <typename DT>
    TritonOutput<DT> fromServer() const;

    //const accessors
    const ShapeView& shape() const { return shape_; }
    int64_t byteSize() const { return byteSize_; }
    const std::string& dname() const { return dname_; }
    unsigned batchSize() const { return batchSize_; }

    //utilities
    bool variableDims() const { return variableDims_; }
    int64_t sizeDims() const { return productDims_; }
    //default to dims if shape isn't filled
    int64_t sizeShape() const { return variableDims_ ? dimProduct(shape_) : sizeDims(); }

  private:
    friend class TritonClient;

    //private accessors only used by client
    bool setShape(const ShapeType& newShape, bool canThrow);
    bool setShape(unsigned loc, int64_t val, bool canThrow);
    void setBatchSize(unsigned bsize);
    void reset();
    void setResult(std::shared_ptr<Result> result) { result_ = result; }
    IO* data() { return data_.get(); }

    //helpers
    bool anyNeg(const ShapeView& vec) const
    {
      return std::any_of(vec.begin(), vec.end(), [](int64_t i) { return i < 0; });
    }
    int64_t dimProduct(const ShapeView& vec) const
    {
      return std::accumulate(vec.begin(), vec.end(), 1, std::multiplies<int64_t>());
    }
    void createObject(IO** ioptr) const;

    //members
    std::string name_;
    std::shared_ptr<IO> data_;
    const ShapeType dims_;
    bool noBatch_;
    unsigned batchSize_;
    ShapeType fullShape_;
    ShapeView shape_;
    bool variableDims_;
    int64_t productDims_;
    std::string dname_;
    inference::DataType dtype_;
    int64_t byteSize_;
    std::any holder_;
    std::shared_ptr<Result> result_;
  };

  using TritonInputData = TritonData<nic::InferInput>;
  using TritonInputMap = std::unordered_map<std::string, TritonInputData>;
  using TritonOutputData = TritonData<nic::InferRequestedOutput>;
  using TritonOutputMap = std::unordered_map<std::string, TritonOutputData>;

  template <>
  void TritonInputData::reset();
  template <>
  void TritonOutputData::reset();
  template <>
  void TritonInputData::createObject(nic::InferInput** ioptr) const;
  template <>
  void TritonOutputData::createObject(nic::InferRequestedOutput** ioptr) const;

  //explicit template instantiation declarations
  extern template class TritonData<nic::InferInput>;
  extern template class TritonData<nic::InferRequestedOutput>;

}
#endif
