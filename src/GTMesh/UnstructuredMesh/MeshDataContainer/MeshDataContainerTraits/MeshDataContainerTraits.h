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
        char name[8] = "dim ";
        nameContainer() {
            unsigned int dim = MeshDataContainer<T, dims...>::template dimensionAt<Index>();
            static_assert (MeshDataContainer<T, dims...>::template dimensionAt<Index>() < 1000, "limit of the string is dimension 999");
            int i = 4;
            for (unsigned int mod = 100; mod > 0; mod /= 10){
                int digit = dim / mod;
                if (digit > 0) {
                    name[i] = '0' + char(digit);
                    i++;
                }
                dim -= digit * mod;
            }
            name[i + 1] = '\0';
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
