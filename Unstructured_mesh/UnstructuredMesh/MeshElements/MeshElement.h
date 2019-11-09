#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H
#include "../UnstructedMeshDefine.h"
#include "../../Vertex.h"
#include "CellBoundaryConnection.h"
#include "ComputationalySignificantElement.h"
#include "type_traits"
#include "tuple"
#include <array>
#include <stdexcept>
#include <vector>

/**
 * @brief The MeshElementBase class provides the basic property
 * of mesh element classes (its index)
 * @see MeshElement
 */
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

    IndexType getIndex(){
        return elementIndex;
    }

    void setIndex(IndexType index){
        elementIndex = index;
    }
/*
    IndexType getGlobalIndex(){
        return globalElementIndex;
    }

    void setGlobalIndex(IndexType index){
        globalElementIndex = index;
    }
*/
};


template <typename IndexType>
struct Subelement{
    IndexType index = INVALID_INDEX(IndexType);
};


template <typename IndexType, unsigned int Reserve>
class SubelementContainer : public std::array<Subelement<IndexType>, Reserve>{
    unsigned int numberOfElements = 0;

public:
    unsigned int getNumberOfSubElements(){
        return numberOfElements;
    }

    unsigned int size() {
        return numberOfElements;
    }

    unsigned int reserve() {
        return Reserve;
    }

    void addSubelement(IndexType index) {
        if (numberOfElements < Reserve){
            this->at(numberOfElements).index = index;
            numberOfElements++;
        } else {
            throw(std::runtime_error(//"In face element (" + std::to_string(MeshElementBase<IndexType>::GetIndex()) +
                                     ") number of edges overgrew the number of reserved indexes (" + std::to_string(Reserve)
                                     +")."));
        }

    }

    void push_back(IndexType index) {
        addSubelement(index);
    }

    void removeSubelement(unsigned char atIndex){
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

    typename std::array<Subelement<IndexType>, Reserve>::iterator end() {
        return this->begin() + getNumberOfSubElements();
    }

    typename std::array<Subelement<IndexType>, Reserve>::const_iterator cend() const {
        return this->cbegin() + getNumberOfSubElements();
    }
};


template<typename IndexType>
class SubelementContainer<IndexType, 0> : public std::vector<IndexType> {
    void addSubelement(IndexType index) {
        this->push_back(index);
    }

    void removeSubelement(unsigned char atIndex){
        this->erase(atIndex);
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

    SubelementContainer<IndexType, Reserve>& getSubelements(){
        return subelements;
    }

    const SubelementContainer<IndexType, Reserve>& getSubelements() const {
        return subelements;
    }

    MeshElement(IndexType index = INVALID_INDEX(IndexType))
        :MeshElementBase<IndexType>(index), CellBoundaryConnection<IndexType> () {}

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
        setVertexAIndex(vertexAIndex);
        setVertexBIndex(vertexBIndex);
    }



    IndexType getVertexAIndex(){
        return vertexAIndex;
    }

    IndexType getVertexBIndex(){
        return vertexBIndex;
    }

    void setVertexAIndex(IndexType index){
        vertexAIndex = index;
    }


    void setVertexBIndex(IndexType index){
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

    IndexType getBoundaryElementIndex(){
        return boundaryElementIndex;
    }

    void setBoundaryElementIndex(IndexType index){
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

    template<unsigned int ElementDimension>
    using ElementType = MeshElement<Dimension, ElementDimension, IndexType, Real, reserve<ElementDimension>()>;


    template<unsigned int dim>
    std::vector<ElementType<dim>>& getElements(){
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return Refs._MeshElements<Dimension, dim, IndexType, Real, Reserve...>::elements;
    }


    template<unsigned int dim>
    const std::vector<ElementType<dim>>&  getElements() const {
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return Refs._MeshElements<Dimension, dim, IndexType, Real, Reserve...>::elements;
    }

    std::vector<Vertex>& getVertices(){
        return getElements<0>();
    }


    std::vector<Edge>& getEdges(){
        return getElements<1>();
    }

    std::vector<Face>& getFaces(){
        return getElements<Dimension - 1>();
    }

    std::vector<Cell>& getCells(){
        return getElements<Dimension>();
    }

    std::vector<Cell>& getBoundaryCells() {
        return BoundaryCells;
    }

private:
    template<unsigned int Dim, typename Dummy = void>
    struct _ClearMesh {
        static void clear(MeshElements<Dimension, IndexType, Real, Reserve...> &mesh){
            mesh.template getElements<Dim>().clear();
            _ClearMesh<Dim - 1>::clear(mesh);
        }
    };


    template<typename Dummy>
    struct _ClearMesh<0, Dummy> {
        static void clear(MeshElements<Dimension, IndexType, Real, Reserve...> &mesh){
            mesh.template getElements<0>().clear();
        }
    };

public:
    /**
     * @brief clear<HR>
     * Sets the size of all vectors in the mesh to 0 by calling clear on them.
     */
    void clear() {
        _ClearMesh<Dimension>::clear(*this);
    }

    void appendBoundaryCell(IndexType cellIndex, IndexType faceIndex){
        Cell c;
        c.setIndex(cellIndex);
        c.SetBoundaryElementIndex(faceIndex);
        BoundaryCells.push_back(c);
    }

    void setupBoundaryCells(){
        for (Face& face : getFaces()){
            if (face.GetCellLeftIndex == INVALID_INDEX(IndexType)){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.SetCellLeftIndex(cellIndex);
                appendBoundaryCell(cellIndex, face.getIndex());
            }
            if (face.GetCellRightIndex == INVALID_INDEX(IndexType)){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.SetCellRightIndex(cellIndex);
                appendBoundaryCell(cellIndex, face.getIndex());
            }
        }
        BoundaryCells.shrink_to_fit();
    }


    void setupBoundaryCellsCenters() {
        for(Cell& cell : BoundaryCells){
            cell.SetCenter(getFaces().at(cell.getBoundaryElementIndex()).GetCenter());
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
        CellSubelementIterator& operator++ () {actual = parentMesh->getFaces().at(actual).getNextBElem(cellIndex) == firstBElem ? INVALID_INDEX(IndexType) : parentMesh->getFaces().at(actual).getNextBElem(cellIndex); return *this;}
        CellSubelementIterator& operator++ (int) {actual = parentMesh->getFaces().at(actual).getNextBElem(cellIndex) == firstBElem ? INVALID_INDEX(IndexType) : parentMesh->getFaces().at(actual).getNextBElem(cellIndex); return *this;}
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
            return CellSubelementIterator(cellIndex, parentMesh->getCells()[cellIndex].getBoundaryElementIndex(), parentMesh);
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
            elementIndex = meshElement.getIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType getIndex(){
            return elementIndex;
        }

        SubelementContainer<IndexType, reserve<ElementDim>()>& getSubelements(){
            return parentMesh->template getElements<ElementDim>()[elementIndex].getSubelements();
        }
    };


    template<typename Dummy>
    class MeshElementWrap<Dimension, Dummy> {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, MeshElement<Dimension, Dimension, IndexType, Real, 0>& meshElement){
            elementIndex = meshElement.getIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType getIndex(){
            return elementIndex;
        }

        CellSubelements getSubelements(){
            return CellSubelements(parentMesh, elementIndex);
        }
    };




    template<typename Dummy>
    class MeshElementWrap<1, Dummy> {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, MeshElement<Dimension, 1, IndexType, Real, 0>& meshElement){
            elementIndex = meshElement.getIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType getIndex(){
            return elementIndex;
        }

        Edge& getElement() {
            return parentMesh->getEdges()[elementIndex];
        }

    };

    template<typename Dummy>
    class MeshElementWrap<0, Dummy> {
        IndexType elementIndex;
        MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh;
    public:
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, Vertex& meshElement){
            elementIndex = meshElement.getIndex();
            this->parentMesh = parentMesh;
        }

        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, IndexType elementIndex){
            this->elementIndex = elementIndex;
            this->parentMesh = parentMesh;
        }

        IndexType getIndex(){
            return elementIndex;
        }

        Vertex& getElement() {
            return parentMesh->GetVertices()[elementIndex];
        }
    };

    template<unsigned int ElementDim>
    MeshElementWrap<ElementDim> getElement(IndexType elementIndex){
        return MeshElementWrap<ElementDim>(this, elementIndex);
    }

};





#endif // MESH_ELEMENT_H

