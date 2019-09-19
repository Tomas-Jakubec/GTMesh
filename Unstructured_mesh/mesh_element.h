#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H
#include "unstructed_mesh_define.h"
#include "vertex.h"
#include "cellboundaryconnection.h"
#include "computationaly_significant_element.h"
#include "type_traits"
#include "tuple"
#include <array>
#include <stdexcept>
#include <vector>

template<typename IndexType>
class MeshElementBase{
    /**
     * @brief Index
     *
     * Index of the element in the mesh or in mesh component
     */
    IndexType ElementIndex;

    /**
     * @brief LocalElementIndex
     *
     * Global index of element in the mesh component
     */
    //IndexType GobalElementIndex;
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
/*
    IndexType GetGlobalIndex(){
        return GlobalElementIndex;
    }

    void SetGlobalIndex(IndexType index){
        GlobalElementIndex = index;
    }
*/
};


template <typename IndexType>
struct Subelement{
    IndexType index = INVALID_INDEX(IndexType);
    bool isLeft = false;
};


template <typename IndexType, unsigned int Reserve>
class SubelementContainer : public std::array<Subelement<IndexType>, Reserve>{
    unsigned char NumberOfElements = 0;

public:
    unsigned char GetNumberOfSubElements(){
        return NumberOfElements;
    }


    void AddSubelement(IndexType index, bool isLeft) {
        if (NumberOfElements < Reserve){
            this->at(NumberOfElements).index = index;
            this->at(NumberOfElements).isLeft = isLeft;
            NumberOfElements++;
        } else {
            throw(std::runtime_error(//"In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") number of edges overgrew the number of reserved indexes (" + std::to_string(Reserve)
                                     +")."));
        }

    }

    void RemoveSubelement(unsigned char atIndex){
        if (atIndex < NumberOfElements){
            for(unsigned char i = atIndex; i < NumberOfElements - 1; i++){
                this->at(i) = this->at(i+1);
            }
            this->at(NumberOfElements) = {INVALID_INDEX(IndexType), false};
            NumberOfElements--;
        } else {
            throw(std::runtime_error(//"In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") removing index " + std::to_string(atIndex)
                                     +" is greather than number of subelements " + std::to_string(NumberOfElements)+ "."));
        }
    }

    typename std::array<Subelement<IndexType>, Reserve>::iterator end(){
        return this->begin() + GetNumberOfSubElements();
    }

    typename std::array<Subelement<IndexType>, Reserve>::const_iterator cend(){
        return this->cbegin() + GetNumberOfSubElements();
    }
};


struct emptyStruct{};

struct emptyStruct2{};





template <unsigned int MeshDim, unsigned int ElementDim, typename IndexType, typename Real, unsigned int Reserve = 0>
class MeshElement : public MeshElementBase<IndexType>,
                    public std::conditional<ElementDim == MeshDim - 1,CellBoundaryConnection<IndexType>, emptyStruct>::type,
                    public std::conditional<ElementDim == MeshDim - 1,ComputationallySignificantElement<MeshDim, Real>, emptyStruct2>::type{
    SubelementContainer<IndexType, Reserve> Subelements;
public:

    SubelementContainer<IndexType, Reserve>& GetSubelements(){
        return Subelements;
    }

    const SubelementContainer<IndexType, Reserve>& GetSubelements() const {
        return Subelements;
    }

    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index), CellBoundaryConnection<IndexType> () {
        Subelements.fill({INVALID_INDEX(IndexType), false});
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







template <unsigned int MeshDim, typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, 0, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>, public Vertex<MeshDim, Real>{
public:
    MeshElement(IndexType index = INVALID_INDEX(IndexType), Vertex<MeshDim, Real> v = {})
        :MeshElementBase<IndexType>(index), Vertex<MeshDim, Real>(v) {

    }
    MeshElement<MeshDim, 0, IndexType, Real>& operator =(std::initializer_list<Real> l){
        this->Vertex<MeshDim,double>::operator=(l);
        return *this;
    }
};





template <unsigned int MeshDim,typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, 1, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>,
          public std::conditional<MeshDim == 2,CellBoundaryConnection<IndexType>, emptyStruct>::type,
          public std::conditional<MeshDim == 2,ComputationallySignificantElement<MeshDim, Real>, emptyStruct2>::type{
public:
    IndexType VertexA;
    IndexType VertexB;
public:

    MeshElement(IndexType index = INVALID_INDEX(IndexType),
                IndexType vertexAIndex = INVALID_INDEX(IndexType),
                IndexType vertexBIndex = INVALID_INDEX(IndexType))
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

};













template <unsigned int MeshDim,typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, MeshDim, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>,
          public ComputationallySignificantElement<MeshDim, Real>{

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

    template<unsigned int dim>
    static unsigned int constexpr reserve(){
        if constexpr (dim == Dimension || dim == 1 || dim == 0){
            return 0;
        } else {
            return std::get<(dim == Dimension || dim == 1 || dim == 0) ? 0 : Dimension - dim - 1>(std::array<unsigned int, sizeof... (Reserve)>{Reserve...});
        }
    }




public:
    _MeshElements<Dimension, Dimension, IndexType, Real, Reserve...> Refs;
    std::vector<MeshElement<Dimension, Dimension, IndexType, Real, 0>> BoundaryCells;

    using Vertex = MeshElement<Dimension, 0, IndexType, Real, 0>;
    using Edge = MeshElement<Dimension, 1, IndexType, Real, 0>;
    using Face = MeshElement<Dimension, Dimension - 1, IndexType, Real, reserve<Dimension - 1>()>;
    using Cell = MeshElement<Dimension, Dimension, IndexType, Real, 0>;

    template <unsigned int dim>
    struct ElemType
    {
        using type = MeshElement<Dimension, dim, IndexType, Real, reserve<dim>()>;
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


    struct CellSubelementIterator: public std::iterator<std::forward_iterator_tag, IndexType>
    {

        IndexType Actual;
        IndexType FirstBElem;
        IndexType Cell;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        CellSubelementIterator(IndexType ci, IndexType act, MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh):Cell(ci){
            FirstBElem = act;
            Actual = FirstBElem;
            this->parentMesh = parentMesh;
        }
        CellSubelementIterator& operator++ () {Actual = parentMesh->GetFaces().at(Actual).GetNextBElem(Cell) == FirstBElem ? INVALID_INDEX(IndexType) : parentMesh->GetFaces().at(Actual).GetNextBElem(Cell); return *this;}
        CellSubelementIterator& operator++ (int) {Actual = parentMesh->GetFaces().at(Actual).GetNextBElem(Cell) == FirstBElem ? INVALID_INDEX(IndexType) : parentMesh->GetFaces().at(Actual).GetNextBElem(Cell); return *this;}
        IndexType operator* (){return Actual;}
        bool operator== (CellSubelementIterator& it) {return Actual == it.Actual;}
        bool operator!= (CellSubelementIterator& it) {return Actual != it.Actual;}
    };

    class CellSubelements {
        IndexType cellIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        CellSubelements(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType cellIndex){
            this->cellIndex = cellIndex;
            this->parentMesh = parentMesh;
        }

        CellSubelementIterator begin() {
            return CellSubelementIterator(cellIndex, parentMesh->GetCells()[cellIndex].GetBoundaryElementIndex(), parentMesh);
        }

        CellSubelementIterator end() {
            return CellSubelementIterator(cellIndex, INVALID_INDEX(IndexType), parentMesh);
        }
    };

public:
    template<unsigned int ElementDim, typename Dummy = void>
    class MeshElementWrap {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, MeshElement<Dimension, ElementDim, IndexType, Real, reserve<ElementDim>()>& meshElement){
            elementIndex = meshElement.GetIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType GetIndex(){
            return elementIndex;
        }

        SubelementContainer<IndexType, reserve<ElementDim>()>& GetSubelements(){
            return parentMesh->template GetElements<ElementDim>()[elementIndex].GetSubelements();
        }
    };


    template<typename Dummy>
    class MeshElementWrap<Dimension, Dummy> {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, MeshElement<Dimension, Dimension, IndexType, Real, 0>& meshElement){
            elementIndex = meshElement.GetIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType GetIndex(){
            return elementIndex;
        }

        CellSubelements GetSubelements(){
            return CellSubelements(parentMesh, elementIndex);
        }
    };




    template<typename Dummy>
    class MeshElementWrap<1, Dummy> {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, MeshElement<Dimension, 1, IndexType, Real, 0>& meshElement){
            elementIndex = meshElement.GetIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType GetIndex(){
            return elementIndex;
        }

        Edge& GetElement() {
            return parentMesh->GetEdges()[elementIndex];
        }

    };

    template<typename Dummy>
    class MeshElementWrap<0, Dummy> {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, Vertex& meshElement){
            elementIndex = meshElement.GetIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType GetIndex(){
            return elementIndex;
        }

        Vertex& GetElement() {
            return parentMesh->GetVertices()[elementIndex];
        }
    };

    template<unsigned int ElementDim>
    MeshElementWrap<ElementDim> GetElement(IndexType elementIndex){
        return MeshElementWrap<ElementDim>(this, elementIndex);
    }

};













#endif // MESH_ELEMENT_H

