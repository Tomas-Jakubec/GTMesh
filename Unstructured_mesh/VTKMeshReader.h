#ifndef VTKMESHREADER_H
#define VTKMESHREADER_H

#include "MeshReader.h"
#include "MeshDataContainer.h"
#include <map>

template<unsigned int MeshDimension>
class VTKMeshReader : public MeshReader<MeshDimension>{

};



template<>
class VTKMeshReader<2> : public MeshReader<2>{
    using reader = MeshReader<2>;
    std::map<int, typename reader::ElementType> TypeConversionTable{
        {3, reader::ElementType::LINE},
        {5, reader::ElementType::TRIANGLE},
        {8, reader::ElementType::QUAD},
        {9, reader::ElementType::QUAD},
        {7, reader::ElementType::POLYGON},
    };

    // file indexing
    //

    //
    //MeshDataContainer<IndexType>
};


#endif // VTKMESHREADER_H
