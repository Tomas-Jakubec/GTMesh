#ifndef MESHREADER_H
#define MESHREADER_H
#include "MeshNativeType.h"

template<unsigned int MeshDimension, typename IndexType, typename Real>
class MeshReader{

};


template <typename IndexType, typename Real>
class MeshReader<2, IndexType, Real> {
public:
    using type = MeshNativeType<2>;


};
#endif // MESHREADER_H
