#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

template<typename T>
using  MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using  VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;


/*
 *
 *  Convert Eigen::Tensor --> Eigen::Matrix
 *
 */


template<typename Derived, typename Device = Eigen::DefaultDevice>
class selfCleaningEvaluator {
private:
    using Evaluator = Eigen::TensorEvaluator<const Eigen::TensorForcedEvalOp<const Derived>, Device>;
    Evaluator m_eval;

public:
    selfCleaningEvaluator(const Evaluator& eval) : m_eval(eval) {}
    selfCleaningEvaluator(const Eigen::TensorBase<Derived, Eigen::ReadOnlyAccessors>& expr, const Device& device = Device())
        : m_eval(Evaluator(expr.eval(), device)) {
        m_eval.evalSubExprsIfNeeded(nullptr);
    }

    ~selfCleaningEvaluator() {
        // This whole point of this object is to call cleanup automatically on destruct.
        // If there are pending operations to evaluate, m_eval will allocate a buffer to hold a result,
        // which needs to be deallocated.
        m_eval.cleanup();
    }

    constexpr auto rank() {
        using DimType = typename decltype(m_eval)::Dimensions::Base;
        return DimType{}.size(); // Because Derived::Dimensions is sometimes wrong
    }
    const Evaluator* operator->() const { return &m_eval; }
    Evaluator* operator->() { return &m_eval; }
    constexpr auto   map() {
        // We inspect m_eval to get the type, rank and dimensions, because it has the resulting tensor,
        // whereas Derived is the tensor type _before_ whatever operation is pending (if any).
        using DimType = typename decltype(m_eval)::Dimensions::Base;
        constexpr auto rank = DimType{}.size();
        using Scalar = typename Eigen::internal::remove_const<typename decltype(m_eval)::Scalar>::type;
        return Eigen::TensorMap<Eigen::Tensor<Scalar, rank>>(m_eval.data(), m_eval.dimensions());
    }
};

// Evaluates expressions if needed
template<typename T, typename Device = Eigen::DefaultDevice>
auto asEval(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr,  // An Eigen::TensorBase object (Tensor, TensorMap, TensorExpr... )
    const Device& device = Device()) {                           // Override to evaluate on another device, e.g. thread pool or gpu.
    return selfCleaningEvaluator(expr, device);
}

// Converts any Eigen::Tensor (or expression) to an Eigen::Matrix with shape rows/cols
template<typename T, typename sizeType, typename Device = Eigen::DefaultDevice>
auto MatrixCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr, const sizeType rows, const sizeType cols, const Device& device = Device()) {
    auto tensor = asEval(expr, device);
    auto tensorMap = tensor.map();
    using Scalar = typename decltype(tensorMap)::Scalar;
    return static_cast<MatrixType<Scalar>>(Eigen::Map<const MatrixType<Scalar>>(tensorMap.data(), rows, cols));
}

// Converts any Eigen::Tensor (or expression) to an Eigen::Matrix with shape rows/cols
template<typename T, typename Device = Eigen::DefaultDevice>
auto VectorCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr, const Device& device = Device()) {
    auto tensor = asEval(expr, device);
    auto tensorMap = tensor.map();
    auto size = Eigen::internal::array_prod(tensorMap.dimensions());
    using Scalar = typename decltype(tensorMap)::Scalar;
    return static_cast<VectorType<Scalar>>(Eigen::Map<const VectorType<Scalar>>(tensorMap.data(), size));
}

// View an existing Eigen::Tensor as an Eigen::Map<Eigen::Matrix>
template<typename Scalar, auto rank, typename sizeType>
auto MatrixMap(const Eigen::Tensor<Scalar, rank>& tensor, const sizeType rows, const sizeType cols) {
    return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
}

// View an existing Eigen::Tensor of rank 2 as an Eigen::Map<Eigen::Matrix>
// Rows/Cols are determined from the matrix
template<typename Scalar>
auto MatrixMap(const Eigen::Tensor<Scalar, 2>& tensor) {
    return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
}

// View an existing Eigen::Tensor of rank 1 as an Eigen::Map<Eigen::Vector>
// Rows is the same as the size of the tensor. 
template<typename Scalar, auto rank>
auto VectorMap(const Eigen::Tensor<Scalar, rank>& tensor) {
    return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
}
