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
     * @brief elementIndex<HR>
     * Index of the element in the mesh or in mesh component
     */
    IndexType elementIndex;

    /**
     * @brief globalElementIndex
     *
     * Global index of element in the mesh component
     */
    //IndexType gobalElementIndex;
public:

    MeshElementBase(IndexType index = INVALID_INDEX(IndexType)){
        elementIndex = index;
    }

    IndexType GetIndex(){
        return elementIndex;
    }

    void SetIndex(IndexType index){
        elementIndex = index;
    }
/*
    IndexType GetGlobalIndex(){
        return globalElementIndex;
    }

    void SetGlobalIndex(IndexType index){
        globalElementIndex = index;
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
    unsigned char numberOfElements = 0;

public:
    unsigned char GetNumberOfSubElements(){
        return numberOfElements;
    }


    void AddSubelement(IndexType index, bool isLeft) {
        if (numberOfElements < Reserve){
            this->at(numberOfElements).index = index;
            this->at(numberOfElements).isLeft = isLeft;
            numberOfElements++;
        } else {
            throw(std::runtime_error(//"In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") number of edges overgrew the number of reserved indexes (" + std::to_string(Reserve)
                                     +")."));
        }

    }

    void RemoveSubelement(unsigned char atIndex){
        if (atIndex < numberOfElements){
            for(unsigned char i = atIndex; i < numberOfElements - 1; i++){
                this->at(i) = this->at(i+1);
            }
            this->at(numberOfElements) = {INVALID_INDEX(IndexType), false};
            numberOfElements--;
        } else {
            throw(std::runtime_error(//"In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") removing index " + std::to_string(atIndex)
                                     +" is greather than number of subelements " + std::to_string(numberOfElements)+ "."));
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
    SubelementContainer<IndexType, Reserve> subelements;
public:

    SubelementContainer<IndexType, Reserve>& GetSubelements(){
        return subelements;
    }

    const SubelementContainer<IndexType, Reserve>& GetSubelements() const {
        return subelements;
    }

    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index), CellBoundaryConnection<IndexType> () {
        subelements.fill({INVALID_INDEX(IndexType), false});
    }

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
    IndexType vertexAIndex;
    IndexType vertexBIndex;
public:

    MeshElement(IndexType index = INVALID_INDEX(IndexType),
                IndexType vertexAIndex = INVALID_INDEX(IndexType),
                IndexType vertexBIndex = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index) {
        SetVertexAIndex(vertexAIndex);
        SetVertexBIndex(vertexBIndex);
    }



    IndexType GetVertexAIndex(){
        return vertexAIndex;
    }

    IndexType GetVertexBIndex(){
        return vertexBIndex;
    }

    void SetVertexAIndex(IndexType index){
        vertexAIndex = index;
    }


    void SetVertexBIndex(IndexType index){
        vertexBIndex = index;
    }

};













template <unsigned int MeshDim,typename IndexType, typename Real, unsigned int Reserve>
class MeshElement<MeshDim, MeshDim, IndexType, Real, Reserve>
        : public MeshElementBase<IndexType>,
          public ComputationallySignificantElement<MeshDim, Real>{

    IndexType boundaryElementIndex;
public:
//GetBoundaryElement
    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index) {

    }

    IndexType GetBoundaryElementIndex(){
        return boundaryElementIndex;
    }

    void SetBoundaryElementIndex(IndexType index){
        boundaryElementIndex = index;
    }

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




    _MeshElements<Dimension, Dimension, IndexType, Real, Reserve...> Refs;
    std::vector<MeshElement<Dimension, Dimension, IndexType, Real, 0>> BoundaryCells;
public:

    using Vertex = MeshElement<Dimension, 0, IndexType, Real, 0>;
    using Edge = MeshElement<Dimension, 1, IndexType, Real, 0>;
    using Face = MeshElement<Dimension, Dimension - 1, IndexType, Real, reserve<Dimension - 1>()>;
    using Cell = MeshElement<Dimension, Dimension, IndexType, Real, 0>;

    template <unsigned int dim>
    struct ElemType //lze nahradit šablonovým using
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

    std::vector<Cell>& GetBoundaryCells() {
        return BoundaryCells;
    }

    void AppendBoundaryCell(IndexType cellIndex, IndexType faceIndex){
        Cell c;
        c.SetIndex(cellIndex);
        c.SetBoundaryElementIndex(faceIndex);
        BoundaryCells.push_back(c);
    }

    void SetupBoundaryCells(){
        for (Face& face : GetFaces()){
            if (face.GetCellLeftIndex == INVALID_INDEX(IndexType)){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.SetCellLeftIndex(cellIndex);
                AppendBoundaryCell(cellIndex, face.GetIndex());
            }
            if (face.GetCellRightIndex == INVALID_INDEX(IndexType)){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.SetCellRightIndex(cellIndex);
                AppendBoundaryCell(cellIndex, face.GetIndex());
            }
        }
        BoundaryCells.shrink_to_fit();
    }


    void SetupBoundaryCellsCenters() {
        for(Cell& cell : BoundaryCells){
            cell.SetCenter(GetFaces().at(cell.GetBoundaryElementIndex()).GetCenter());
        }
    }

    struct CellSubelementIterator: public std::iterator<std::forward_iterator_tag, IndexType>
    {

        IndexType actual;
        IndexType firstBElem;
        IndexType cellIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        CellSubelementIterator(IndexType ci, IndexType act, MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh):cellIndex(ci){
            firstBElem = act;
            actual = firstBElem;
            this->parentMesh = parentMesh;
        }
        CellSubelementIterator& operator++ () {actual = parentMesh->GetFaces().at(actual).GetNextBElem(cellIndex) == firstBElem ? INVALID_INDEX(IndexType) : parentMesh->GetFaces().at(actual).GetNextBElem(cellIndex); return *this;}
        CellSubelementIterator& operator++ (int) {actual = parentMesh->GetFaces().at(actual).GetNextBElem(cellIndex) == firstBElem ? INVALID_INDEX(IndexType) : parentMesh->GetFaces().at(actual).GetNextBElem(cellIndex); return *this;}
        IndexType operator* (){return actual;}
        bool operator== (CellSubelementIterator& it) {return actual == it.actual;}
        bool operator!= (CellSubelementIterator& it) {return actual != it.actual;}
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

