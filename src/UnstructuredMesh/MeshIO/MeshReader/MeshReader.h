#ifndef MESHREADER_H
#define MESHREADER_H
#include "../MeshNativeType.h"
#include "../../MeshDataContainer/MeshDataContainer.h"

template<unsigned int MeshDimension>
class MeshReader{
public:
    using type = MeshNativeType<MeshDimension>;

    virtual
    MeshDataContainer<typename type::ElementType, MeshDimension>
    getCellTypes() const = 0;
};


#endif // MESHREADER_H
