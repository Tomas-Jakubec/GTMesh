#ifndef UNSTRUCTED_MESH_DEFINE_H
#define UNSTRUCTED_MESH_DEFINE_H

#include <limits>
#define INVALID_INDEX(IndexType) (std::numeric_limits<IndexType>::max())
#define BOUNDARY_INDEX(IndexType) (static_cast<IndexType>(1) << (std::numeric_limits<IndexType>::digits - 1))
#define EXTRACTING_INDEX(IndexType) (static_cast<IndexType>(~BOUNDARY_INDEX(IndexType)))


template <typename IndexType>
bool isInvalidIndex(const IndexType& index){
    return index == INVALID_INDEX(IndexType);
}


template <typename IndexType>
bool isBoundaryIndex(const IndexType& index){
    return (BOUNDARY_INDEX(IndexType) & index) == BOUNDARY_INDEX(IndexType);
}


template <typename IndexType>
IndexType makeBoundaryIndex(const IndexType& index){
    return (EXTRACTING_INDEX(IndexType) | index);
}


template <typename IndexType>
IndexType extractBoundaryIndex(const IndexType& index){
    return (EXTRACTING_INDEX(IndexType) & index);
}

#endif // UNSTRUCTED_MESH_DEFINE_H
