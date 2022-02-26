#ifndef ELEMENTINDEX_HPP
#define ELEMENTINDEX_HPP

#include <cstddef>

/**
 * @brief Small integer wrapper to bound a element dimension with an index.
 */
template<std::size_t ElementDimension, typename IndexType = std::size_t>
struct ElementIndex {
    IndexType index;

    ElementIndex(const IndexType& oIndex){
        index = oIndex;
    }

    ElementIndex(const ElementIndex&) = default;

    ElementIndex operator=(const IndexType& rhs){
        index = rhs;
    }

    operator IndexType() {
        return index;
    }
};

#endif // ELEMENTINDEX_HPP
