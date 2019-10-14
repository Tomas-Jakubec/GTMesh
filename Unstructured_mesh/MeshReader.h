#ifndef MESHREADER_H
#define MESHREADER_H
#include "MeshNativeType.h"

template<unsigned int MeshDimension>
class MeshReader{
public:
    using type = MeshNativeType<MeshDimension>;
};


#endif // MESHREADER_H
