#ifndef MESHWRITER_H
#define MESHWRITER_H
#include "MeshNativeType.h"

template<unsigned int MeshDimension, typename IndexType, typename Real>
class MeshWriter{

};


template <typename IndexType, typename Real>
class MeshWriter<2, IndexType, Real> {
public:
    using type = MeshNativeType<2>;

};

template <typename IndexType, typename Real>
class MeshWriter<3, IndexType, Real> {
public:
    using type = MeshNativeType<3>;

};
#endif // MESHWRITER_H
