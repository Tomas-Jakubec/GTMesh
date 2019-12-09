#ifndef MESHELEMENTS_H
#define MESHELEMENTS_H

#include "MeshElement.h"

template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
struct MeshElements{
private:

    template<unsigned int dim, typename Void = void>
    struct _Reserve{

        static unsigned int constexpr value = std::get<Dimension - dim - 1>(std::array<unsigned int, sizeof... (Reserve)>{Reserve...});

    };

    template<unsigned int dim>
    struct _Reserve<dim, typename std::enable_if<dim == Dimension || dim == 1 || dim == 0 || (Dimension - dim > sizeof...(Reserve))>::type>{

        static unsigned int constexpr value = 0;
    };

public:

    template<unsigned int dim>
    static unsigned int constexpr reserve() {
        return _Reserve<dim>::value;
    }

    using Vertex = MeshElement<Dimension, 0, IndexType, Real, 0>;
    using Edge = MeshElement<Dimension, 1, IndexType, Real, 0>;
    using Face = MeshElement<Dimension, Dimension - 1, IndexType, Real, _Reserve<Dimension - 1>::value>;
    using Cell = MeshElement<Dimension, Dimension, IndexType, Real, 0>;

    template<unsigned int ElementDimension>
    using ElementType = MeshElement<Dimension, ElementDimension, IndexType, Real, _Reserve<ElementDimension>::value>;

private:
    template <unsigned int ElemDim = Dimension, typename Dummy = void>
    struct _MeshElements : public _MeshElements<ElemDim - 1, Dummy>{
        std::vector<typename MeshElements<Dimension, IndexType, Real, Reserve...>:: template ElementType<ElemDim>> elements;
    };

    template <typename Dummy>
    struct _MeshElements<0, Dummy>{
        std::vector<typename MeshElements<Dimension, IndexType, Real, Reserve...>:: template ElementType<0>> elements;
    };




private:
    _MeshElements<Dimension> innerElements;
    std::vector<Cell> BoundaryCells;
    //_MeshElements<Dimension> boundaryElements;


public:
    template<unsigned int dim>
    std::vector<ElementType<dim>>& getElements(){
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return innerElements._MeshElements<dim>::elements;
    }


    template<unsigned int dim>
    const std::vector<ElementType<dim>>&  getElements() const {
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return innerElements._MeshElements<dim>::elements;
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
/*
 * Constant version of getters
 */
    const std::vector<Vertex>& getVertices() const {
        return getElements<0>();
    }


    const std::vector<Edge>& getEdges() const {
        return getElements<1>();
    }

    const std::vector<Face>& getFaces() const {
        return getElements<Dimension - 1>();
    }

    const std::vector<Cell>& getCells() const {
        return getElements<Dimension>();
    }

    const std::vector<Cell>& getBoundaryCells() const {
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
        c.setBoundaryElementIndex(faceIndex);
        BoundaryCells.push_back(c);
    }

    void setupBoundaryCells(){
        for (Face& face : getFaces()){
            if (face.getCellLeftIndex() == INVALID_INDEX(IndexType)){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.setCellLeftIndex(cellIndex);
                appendBoundaryCell(cellIndex, face.getIndex());
            }
            if (face.getCellRightIndex() == INVALID_INDEX(IndexType)){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.setCellRightIndex(cellIndex);
                appendBoundaryCell(cellIndex, face.getIndex());
            }
        }
        BoundaryCells.shrink_to_fit();
    }


    void setupBoundaryCellsCenters() {
        for(Cell& cell : BoundaryCells){
            cell.setCenter(getFaces().at(cell.getBoundaryElementIndex()).getCenter());
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
        MeshElementWrap(MeshElements<Dimension, IndexType, Real, Reserve...>* parentMesh, MeshElement<Dimension, ElementDim, IndexType, Real, _Reserve<Dimension - 1>::value>& meshElement){
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

        SubelementContainer<IndexType, _Reserve<Dimension - 1>::value>& getSubelements(){
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


#endif // MESHELEMENTS_H
