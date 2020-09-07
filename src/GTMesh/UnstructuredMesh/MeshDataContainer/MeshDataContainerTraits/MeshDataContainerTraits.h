#ifndef MESHDATACONTAINERTRAITS_H
#define MESHDATACONTAINERTRAITS_H

#include "../MeshDataContainer.h"
#include <GTMesh/Traits/Traits.h>

template <typename T, unsigned int ...dims>
class Traits<MeshDataContainer<T, dims...>> {
public:
    template< unsigned int Index = 0,
              bool = (Index == sizeof... (dims)) > // stopping condition
    struct nameContainer : nameContainer<Index + 1>{
        char name[9] = {};
        nameContainer() {
            sprintf(name, "dim: %d", MeshDataContainer<T, dims...>::template dimensionAt<Index>());
        }
    };

    template<unsigned int Index>
    struct nameContainer<Index, true> {};
    nameContainer<> names;
public:

    template<unsigned int Index>
    using type = typename MeshDataContainer<T, dims...>::template DataContainerType<Index>;

    static constexpr unsigned int size() {
        return sizeof... (dims);
    }

    template <unsigned int Index>
    auto& getAttr(MeshDataContainer<T, dims...>& t) const {
        return t.template getDataByPos<Index>();
    }

    template <unsigned int Index>
    const auto& getValue(const MeshDataContainer<T, dims...>& t) const {
        return t.template getDataByPos<Index>();
    }


    template <unsigned int Index>
    void setValue(
            MeshDataContainer<T, dims...>& t,
            const type<Index>& val) const {
        t.template getDataByPos<Index>() = val;
    }

    template <unsigned int Index>
    const char* getName() const {
        return names.nameContainer<Index, false>::name;
    }
};


#endif // MESHDATACONTAINERTRAITS_H
