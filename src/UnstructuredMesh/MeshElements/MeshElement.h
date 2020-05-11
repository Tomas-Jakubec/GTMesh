#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H
#include "../UnstructuredMeshDefine.h"
#include "../../NumericStaticArray/Vertex.h"
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

    IndexType getIndex() const {
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





template <typename IndexType, unsigned int Reserve>
class SubelementContainer : public std::array<IndexType, Reserve>{
    unsigned int numberOfElements = 0;

public:
    unsigned int getNumberOfSubElements() const {
        return numberOfElements;
    }

    unsigned int size() const {
        return numberOfElements;
    }

    unsigned int reserve() const {
        return Reserve;
    }

    void addSubelement(IndexType index) {
        if (numberOfElements < Reserve){
            this->at(numberOfElements) = index;
            numberOfElements++;
        } else {
            throw(
                std::runtime_error(
                        "number of edges overgrew the number of reserved indexes (" +
                        std::to_string(Reserve) + ")."
                        )
                 );
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
            throw(std::runtime_error(
                        "removing index " + std::to_string(atIndex) +
                        " is greather than number of subelements " + std::to_string(numberOfElements)+ "."
                        )
                 );
        }
    }

    typename std::array<IndexType, Reserve>::iterator end() {
        return this->begin() + getNumberOfSubElements();
    }

    typename std::array<IndexType, Reserve>::const_iterator end() const {
        return this->cbegin() + getNumberOfSubElements();
    }

    typename std::array<IndexType, Reserve>::const_iterator cend() const {
        return this->cbegin() + getNumberOfSubElements();
    }
};


template<typename IndexType>
class SubelementContainer<IndexType, 0> : public std::vector<IndexType> {
public:
    IndexType getNumberOfSubElements() const {
        return this->size();
    }

    void addSubelement(IndexType index) {
        this->push_back(IndexType{index});
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



    IndexType getVertexAIndex() const {
        return vertexAIndex;
    }

    IndexType getVertexBIndex() const {
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

    IndexType getBoundaryElementIndex() const {
        return boundaryElementIndex;
    }

    void setBoundaryElementIndex(IndexType index){
        boundaryElementIndex = index;
    }

};




#endif // MESH_ELEMENT_H

