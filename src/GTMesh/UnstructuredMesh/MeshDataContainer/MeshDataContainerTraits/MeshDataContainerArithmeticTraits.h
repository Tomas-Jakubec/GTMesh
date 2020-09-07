#ifndef MESHDATACONTAINERARITHMETICTRAITS_H
#define MESHDATACONTAINERARITHMETICTRAITS_H
#include "MeshDataContainerTraits.h"

template <typename T, unsigned int ...dims>
class DefaultArithmeticTraits<MeshDataContainer<T, dims...>> {
public:
    using traitsType = Traits<MeshDataContainer<T, dims...>>;

    static traitsType getTraits(){
        return traitsType();
    }

    static constexpr unsigned int size() {
        return traitsType::size();
    }
};

#endif // MESHDATACONTAINERARITHMETICTRAITS_H
