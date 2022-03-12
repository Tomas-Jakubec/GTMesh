#ifndef ELEMENTINDEX_HPP
#define ELEMENTINDEX_HPP

#include <cstddef>

/**
 * @brief Small integer wrapper to bound a element dimension with an index.
 */
template<unsigned int ElementDimension, typename IndexType = unsigned int>
struct ElementIndex {
    IndexType index;

    ElementIndex(const ElementIndex&) = default;
    ElementIndex& operator=(const ElementIndex& rhs){
        index = rhs.index;
        return *this;
    }

    ElementIndex(const IndexType& oIndex){
        index = oIndex;
    }

    ElementIndex& operator=(const IndexType& rhs){
        index = rhs;
        return *this;
    }

    bool operator==(const ElementIndex& rhs) const {
        return index == rhs.index;
    }

    operator IndexType() {
        return index;
    }

    static constexpr unsigned int getDimension(){
        return ElementDimension;
    }
};

#endif // ELEMENTINDEX_HPP
