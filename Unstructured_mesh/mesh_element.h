#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H
#include "unstructed_mesh_define.h"
#include "vertex.h"
#include "cellboundaryconnection.h"
#include "type_traits"
template<typename indexType>
class MeshElementBase{
    /**
     * @brief Index
     *
     * Index of the element in the mesh
     */
    indexType ElementIndex;
public:

    MeshElementBase(indexType index = INVALID_INDEX(indexType)){
        ElementIndex = index;
    }

    indexType GetIndex(){
        return ElementIndex;
    }

    void SetIndex(indexType index){
        ElementIndex = index;
    }

};

class emptyClass{

};






template <unsigned int MeshDim, unsigned int ElementDim, typename indexType, typename Real>
class MeshElement : public MeshElementBase<indexType>{

};






template <unsigned int MeshDim, typename indexType, typename Real>
class MeshElement<MeshDim, 0, indexType, Real> : public MeshElementBase<indexType>, public Vertex<MeshDim, Real>{
public:
    MeshElement(indexType index = INVALID_INDEX(indexType))
        :MeshElementBase<indexType>(index), Vertex<MeshDim, Real>() {

    }
    MeshElement<MeshDim, 0, indexType, Real>& operator =(std::initializer_list<Real> l){
        this->Vertex<3,double>::operator=(l);
        return *this;
    }
};





template <unsigned int MeshDim,typename indexType, typename Real>
class MeshElement<MeshDim, 1, indexType, Real>
        : public MeshElementBase<indexType>,
        public std::conditional<MeshDim == 2,CellBoundaryConnection<indexType>, emptyClass>::type{

    indexType VertexA;
    indexType VertexB;
public:
  //TODO doplnit metody pro práci hranami
    MeshElement(indexType index = INVALID_INDEX(indexType))
        :MeshElementBase<indexType>(index), std::conditional<MeshDim == 2,CellBoundaryConnection<indexType>, emptyClass>::type() {

    }

};


template <typename indexType, typename Real>
class MeshElement<3, 2, indexType, Real>
        : public MeshElementBase<indexType>,
        public CellBoundaryConnection<indexType>{

    indexType subElements[6];
    unsigned char numberOfElements = 0;
public:
  //TODO doplnit metody pro práci hranami
    MeshElement(indexType index = INVALID_INDEX(indexType))
        :MeshElementBase<indexType>(index), CellBoundaryConnection<indexType> () {

    }

};


template <unsigned int MeshDim,typename indexType, typename Real>
class MeshElement<MeshDim, MeshDim, indexType, Real>
        : public MeshElementBase<indexType>{

    indexType BoundaryElement;
    Vertex<MeshDim, Real> Center;
public:
//GetBoundaryElement
    MeshElement(indexType index = INVALID_INDEX(indexType))
        :MeshElementBase<indexType>(index) {

    }
    const Vertex<MeshDim, Real>& CalculateCenter();
};




#endif // MESH_ELEMENT_H
