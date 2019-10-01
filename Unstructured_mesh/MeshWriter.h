#ifndef MESHWRITER_H
#define MESHWRITER_H
#include "MeshNativeType.h"

template<unsigned int MeshDimension, typename IndexType, typename Real>
class MeshWriter{

};

enum {
    hmmm
};

template <typename IndexType, typename Real>
class MeshWriter<2, IndexType, Real> {
public:
    using type = MeshNativeType<2>;

};
#endif // MESHWRITER_H
