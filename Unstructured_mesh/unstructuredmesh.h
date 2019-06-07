#ifndef UNSTRUCTUREDMESH_H
#define UNSTRUCTUREDMESH_H
#include "mesh_element.h"
#include "vector"

template <unsigned int Dimension,unsigned int ElemDim, typename indexType, typename Real>
struct MeshElements : public MeshElements<Dimension, ElemDim - 1, indexType, Real>{
    std::vector<MeshElement<Dimension, ElemDim, indexType, Real>> elements;
};

template <unsigned int Dimension, typename indexType, typename Real>
struct MeshElements<Dimension, 0, indexType, Real>{
    std::vector<MeshElement<Dimension, 0, indexType, Real>> elements;
};


template <unsigned int Dimension, typename indexType, typename Real>
class UnstructuredMesh{
    MeshElements<Dimension, Dimension, indexType, Real> elements;
public:

    using Cell = MeshElement<Dimension, Dimension, indexType, Real>;

    template<unsigned int dim>
    std::vector<MeshElement<Dimension, dim, indexType, Real>>& GetElements(){
        return elements.MeshElements<Dimension, dim, indexType, Real>::elements;
    }

    std::vector<MeshElement<Dimension, Dimension, indexType, Real>>& GetCells(){
        return GetElements<Dimension>();
    }


    std::vector<MeshElement<Dimension, 0, indexType, Real>>& GetPoints(){
        return GetElements<0>();
    }

    std::vector<MeshElement<Dimension, 0, indexType, Real>>& GetPoints(){
        return GetElements<0>();
    }
};

#endif // UNSTRUCTUREDMESH_H
