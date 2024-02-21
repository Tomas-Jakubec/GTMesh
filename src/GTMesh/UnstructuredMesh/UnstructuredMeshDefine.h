#ifndef UNSTRUCTED_MESH_DEFINE_H
#define UNSTRUCTED_MESH_DEFINE_H

#include <limits>
#define INVALID_INDEX(IndexType) (std::numeric_limits<IndexType>::max())
#define BOUNDARY_INDEX(IndexType) (static_cast<IndexType>(1) << (std::numeric_limits<IndexType>::digits - 1))
#define EXTRACTING_INDEX(IndexType) (static_cast<IndexType>(~BOUNDARY_INDEX(IndexType)))


template <typename IndexType>
constexpr IndexType invalidIndex() {
    return std::numeric_limits<IndexType>::max();
}

template <typename IndexType>
constexpr IndexType boundaryIndexBase() {
    return static_cast<IndexType>(1) << (std::numeric_limits<IndexType>::digits - 1);
}

template <typename IndexType>
constexpr IndexType boundaryExtractIndexBase() {
    return static_cast<IndexType>(~boundaryIndexBase<IndexType>());
}

template <typename IndexType>
bool isInvalidIndex(const IndexType& index){
    return index == invalidIndex<IndexType>();
}


template <typename IndexType>
bool isBoundaryIndex(const IndexType& index){
    return (boundaryIndexBase<IndexType>() & index) == boundaryIndexBase<IndexType>();
}


template <typename IndexType>
IndexType makeBoundaryIndex(const IndexType& index){
    return (boundaryIndexBase<IndexType>() | index);
}


template <typename IndexType>
IndexType extractBoundaryIndex(const IndexType& index){
    return (boundaryExtractIndexBase<IndexType>() & index);
}

#endif // UNSTRUCTED_MESH_DEFINE_H
