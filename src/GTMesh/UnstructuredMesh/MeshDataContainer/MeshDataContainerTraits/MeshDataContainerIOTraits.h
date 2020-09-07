#ifndef MESHDATACONTAINERIOTRAITS_H
#define MESHDATACONTAINERIOTRAITS_H
#include "MeshDataContainerTraits.h"
template <typename T, unsigned int ...dims>
class DefaultIOTraits<MeshDataContainer<T, dims...>> {
public:
    using traitsType = Traits<MeshDataContainer<T, dims...>>;

    static traitsType getTraits(){
        return traitsType();
    }

    static constexpr unsigned int size() {
        return traitsType::size();
    }
};
#endif // MESHDATACONTAINERIOTRAITS_H
