#ifndef UNSTRUCTED_MESH_DEFINE_H
#define UNSTRUCTED_MESH_DEFINE_H

#include <limits>
#define INVALID_INDEX(indexType) (std::numeric_limits<indexType>::max())
#define BOUNDARY_INDEX(indexType) (static_cast<indexType>(1) << (std::numeric_limits<indexType>::digits - 1))
#define EXTRACTING_INDEX(indexType) (static_cast<indexType>(~BOUNDARY_INDEX(indexType)))

#endif // UNSTRUCTED_MESH_DEFINE_H