#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H
#include "unstructed_mesh_define.h"
#include "vertex.h"
#include "cellboundaryconnection.h"
#include "computationaly_significant_element.h"
#include "type_traits"
#include "tuple"
#include <stdexcept>
#include <vector>

template<typename IndexType>
class MeshElementBase{
    /**
     * @brief Index
     *
     * Index of the element in the mesh
     */
    IndexType ElementIndex;
public:

    MeshElementBase(IndexType index = INVALID_INDEX(IndexType)){
        ElementIndex = index;
    }

    IndexType GetIndex(){
        return ElementIndex;
    }

    void SetIndex(IndexType index){
        ElementIndex = index;
    }

};








template <unsigned int MeshDim, unsigned int ElementDim, typename IndexType, typename Real, unsigned int Reserve = 0>
class MeshElement : public MeshElementBase<IndexType>{

};










template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
struct MeshElements{
private:
    template <unsigned int _Dimension,unsigned int ElemDim, typename _IndexType, typename _Real, unsigned int ..._Reserve>
    struct _MeshElements : public _MeshElements<_Dimension, ElemDim - 1, _IndexType, _Real, _Reserve...>{
        std::vector<MeshElement<_Dimension, ElemDim, _IndexType, _Real, std::get<_Dimension - ElemDim - 1>(std::array<unsigned int, sizeof... (Reserve)>{Reserve...})>> elements;
    };

    template <unsigned int _Dimension, typename _IndexType, typename _Real, unsigned int ..._Reserve>
    struct _MeshElements<_Dimension, 0, _IndexType, _Real, _Reserve...>{
        std::vector<MeshElement<Dimension, 0, IndexType, Real, 0>> elements;
    };

    template <unsigned int _Dimension, typename _IndexType, typename _Real, unsigned int ..._Reserve>
    struct _MeshElements<_Dimension, 1, _IndexType, _Real, _Reserve...> : public _MeshElements<_Dimension, 0, _IndexType, _Real, _Reserve...>{
        std::vector<MeshElement<Dimension, 1, IndexType, Real, 0>> elements;
    };


    template <unsigned int _Dimension,typename _IndexType, typename _Real, unsigned int ..._Reserve>
    struct _MeshElements<_Dimension, _Dimension, _IndexType, _Real, _Reserve...> : public _MeshElements<_Dimension, _Dimension - 1, _IndexType, _Real, _Reserve...>{
        std::vector<MeshElement<_Dimension, _Dimension, _IndexType, _Real, 0>> elements;
    };


public:
    _MeshElements<Dimension, Dimension, IndexType, Real, Reserve...> Refs;

    using Vertex = MeshElement<Dimension, 0, IndexType, Real, 0>;
    using Edge = MeshElement<Dimension, 1, IndexType, Real, 0>;
    using Face = MeshElement<Dimension, Dimension - 1, IndexType, Real, Dimension <= 2 ? 0 : std::get<0>(std::array<unsigned int, sizeof... (Reserve)>{Reserve...})>;
    using Cell = MeshElement<Dimension, Dimension, IndexType, Real, 0>;

    template <unsigned int dim>
    struct ElemType
    {
        using type = MeshElement<Dimension, dim, IndexType, Real, (dim == Dimension || dim == 1 || dim == 0) ? 0 : std::get<(dim == Dimension || dim == 1 || dim == 0) ? 0 : Dimension - dim - 1>(std::array<unsigned int, sizeof... (Reserve)>{Reserve...})>;
    };

    template<unsigned int dim>
    std::vector<typename ElemType<dim>::type>& GetElements(){
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return Refs._MeshElements<Dimension, dim, IndexType, Real, Reserve...>::elements;
    }


    template<unsigned int dim>
    const std::vector<typename ElemType<dim>::type>&  GetElements() const {
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return Refs._MeshElements<Dimension, dim, IndexType, Real, Reserve...>::elements;
    }

    std::vector<Vertex>& GetVertices(){
        return GetElements<0>();
    }


    std::vector<Edge>& GetEdges(){
        return GetElements<1>();
    }

    std::vector<Face>& GetFaces(){
        return GetElements<Dimension - 1>();
    }

    std::vector<Cell>& GetCells(){
        return GetElements<Dimension>();
    }

};







template <unsigned int MeshDim, typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, 0, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>, public Vertex<MeshDim, Real>{
public:
    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index), Vertex<MeshDim, Real>() {

    }
    MeshElement<MeshDim, 0, IndexType, Real>& operator =(std::initializer_list<Real> l){
        this->Vertex<MeshDim,double>::operator=(l);
        return *this;
    }
};





template <unsigned int MeshDim,typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, 1, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>{
public:
    IndexType VertexA;
    IndexType VertexB;
public:

    MeshElement(IndexType vertexAIndex = INVALID_INDEX(IndexType),
                IndexType vertexBIndex = INVALID_INDEX(IndexType),
                IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index) {
        SetVertexAIndex(vertexAIndex);
        SetVertexBIndex(vertexBIndex);
    }



    IndexType GetVertexAIndex(){
        return VertexA;
    }

    IndexType GetVertexBIndex(){
        return VertexB;
    }

    void SetVertexAIndex(IndexType index){
        VertexA = index;
    }


    void SetVertexBIndex(IndexType index){
        VertexB = index;
    }
/*
    Real CalculateMeasure(const MeshElements<MeshDim, IndexType, Real, Reserve>& sube){
        auto& vertices = sube.template GetElements<0>();
        return (vertices.elements->at(VertexA) - vertices.elements->at(VertexB)).NormEukleid();
    }

    Vertex<MeshDim, Real> ComputeCenter(const MeshElements<MeshDim, IndexType, Real, Reserve>& sube){
        auto& vertices = sube.template GetElements<0>();
        return (vertices.at(VertexA) + vertices.at(VertexB)) * 0.5;
    }*/
};



template <typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<2, 1, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>,
        public CellBoundaryConnection<IndexType>,
        public ComputationalySignificantElement<2, Real>{
public:
    IndexType VertexA;
    IndexType VertexB;
public:


    MeshElement(IndexType vertexAIndex = INVALID_INDEX(IndexType),
                IndexType vertexBIndex = INVALID_INDEX(IndexType),
                IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index) {
        SetVertexAIndex(vertexAIndex);
        SetVertexBIndex(vertexBIndex);
    }



    IndexType GetVertexAIndex(){
        return VertexA;
    }

    IndexType GetVertexBIndex(){
        return VertexB;
    }

    void SetVertexAIndex(IndexType index){
        VertexA = index;
    }


    void SetVertexBIndex(IndexType index){
        VertexB = index;
    }

    Real CalculateMeasure(const MeshElements<2,IndexType, Real, Reserve>& ref){
        return (ref.template GetVector<0>().at(VertexA) - ref.template GetVector<0>().at(VertexB)).NormEukleid();
    }
/*
    Real CalculateMeasureOverCellDist(const MeshElements<2,IndexType, Real, Reserve>& ref){
        return CalculateMeasure() / (ref.template GetVector<2>().at(this->GetCellLeftIndex()).GetCenter() - ref.template GetVector<2>().at(this->GetCellRightIndex()).GetCenter()).NormEukleid();
    }

    Vertex<2, Real> ComputeCenter(const MeshElements<2,IndexType, Real, Reserve>& ref){
        auto& vertices = ref.template GetElements<0>();
        ComputationalySignificantElement<2, Real>::Center = (vertices.at(VertexA) + vertices.at(VertexB)) * 0.5;
        return ComputationalySignificantElement<2, Real>::Center;
    }
*/
};








template <typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<3, 2, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>,
        public CellBoundaryConnection<IndexType>,
        public ComputationalySignificantElement<3, Real>{

    IndexType subElements[Reserve];
    unsigned char numberOfElements = 0;
public:

    unsigned char GetNumberOfSubElements(){
        return numberOfElements;
    }

    IndexType GetSubelementIndex(unsigned char index){
        return subElements[index];
    }

    void AddSubelementIndex(IndexType index) {
        if (numberOfElements < Reserve){
            subElements[numberOfElements] = index;
            numberOfElements++;
        } else {
            throw(std::runtime_error("In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") number of edges overgrew the number of reserved indexes (" + std::to_string(Reserve)
                                     +")."));
        }

    }

    void RemoveSubelementIndex(unsigned char atIndex){
        if (atIndex < numberOfElements){
            for(unsigned char i = atIndex; i < numberOfElements - 1; i++){
                subElements[i] = subElements[i+1];
            }
            subElements[numberOfElements] = INVALID_INDEX(IndexType);
            numberOfElements--;
        } else {
            throw(std::runtime_error("In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") removing index " + std::to_string(atIndex)
                                     +" is greather than number of subelements " + std::to_string(numberOfElements)+ "."));
        }
    }

    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index), CellBoundaryConnection<IndexType> () {
        for (unsigned char i = 0; i < Reserve; i++) {
            subElements[i] = INVALID_INDEX(IndexType);
        }
    }
/*
    Vertex<3, Real> ComputeCenter(const MeshElements<3, IndexType, Real, Reserve>& ref){

        Vertex<3,Real> tmpCenter = {};
        for(unsigned char i = 0; i < numberOfElements; i++){
            tmpCenter += ref.template GetElements<1>().at(subElements[i]).ComputeCenter(ref);
        }
        ComputationalySignificantElement<3, Real>::Center = tmpCenter * (1.0 / numberOfElements);
    }
    */
};


template <unsigned int MeshDim,typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, MeshDim, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>,
          public ComputationalySignificantElement<MeshDim, Real>{

    IndexType BoundaryElement;
public:
//GetBoundaryElement
    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index) {

    }

    IndexType GetBoundaryElementIndex(){
        return BoundaryElement;
    }

    void SetBoundaryElementIndex(IndexType index){
        BoundaryElement = index;
    }
/*
    Vertex<MeshDim, Real> ComputeCenter(MeshElements<MeshDim, IndexType, Real, Reserve>& ref){
        IndexType boundary = this->BoundaryElement;
        size_t tmp_boundary = boundary;
        Vertex<MeshDim, Real> tmpCenter = {};
        int numberOfBoundaries = 0;
        do {
            numberOfBoundaries++;
            tmpCenter += ref.GetFaces().at(tmp_boundary).ComputeCenter(ref);
            tmp_boundary = ref.GetFaces().at(tmp_boundary).GetNextBElem(this->GetIndex());
        } while(boundary != tmp_boundary);
        ComputationalySignificantElement<MeshDim, Real>::Center = tmpCenter * (1.0 / numberOfBoundaries);
        return ComputationalySignificantElement<MeshDim, Real>::Center;
    }
*/
};




#endif // MESH_ELEMENT_H

