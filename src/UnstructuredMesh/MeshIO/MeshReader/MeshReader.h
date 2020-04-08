#ifndef MESHREADER_H
#define MESHREADER_H
#include "../MeshNativeType.h"
#include "../../MeshDataContainer/MeshDataContainer.h"

/**
 * The base class of MeshReaders
 */
template<unsigned int MeshDimension>
class MeshReader{
public:
    using type = MeshNativeType<MeshDimension>;

    /**
     * @brief Returns the types of the loaded cells.
     */
    virtual
    MeshDataContainer<typename type::ElementType, MeshDimension>
    getCellTypes() const = 0;
};


#endif // MESHREADER_H
