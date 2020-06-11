#ifndef CELLBOUNDARYCONNECTION_H
#define CELLBOUNDARYCONNECTION_H
#include "../UnstructuredMeshDefine.h"
#include "CellConnection.h"
#include <stdexcept>
#include <sstream>

template<typename IndexType>
class CellBoundaryConnection : public CellConnection<IndexType> {
    /**
     * @brief nextBElemWRTCL
     *
     * Index of the next boundary element with respect to the left cell
     */
    IndexType nextBElemWRTCL;
    /**
     * @brief nextBElemWRTCR
     *
     * Index of the next boundary element with respect to right cell
     */
    IndexType nextBElemWRTCR;
public:
    CellBoundaryConnection(IndexType cellLeft = INVALID_INDEX(IndexType),
                           IndexType cellRight = INVALID_INDEX(IndexType),
                           IndexType nextLeft = INVALID_INDEX(IndexType),
                           IndexType nextRight = INVALID_INDEX(IndexType));
    /*
    ** Set atribute methods
    */
    void setNextBElemWRTCR(IndexType nextIndex);

    void setNextBElemWRTCL(IndexType nextIndex);

    bool setNextBElem(IndexType nextBElemIndex, IndexType cellIndex);

    IndexType getNextBElemWRTCL() const;

    IndexType getNextBElemWRTCR() const;

    IndexType getNextBElem(IndexType cellIndex) const;


};

/**
 * Implementation of methods
 */

template<typename IndexType>
CellBoundaryConnection<IndexType>::CellBoundaryConnection(IndexType cellLeft, IndexType cellRight, IndexType nextLeft, IndexType nextRight)
    : CellConnection<IndexType> (cellLeft, cellRight) {
    nextBElemWRTCL = nextLeft;
    nextBElemWRTCR = nextRight;
}

template<typename IndexType>
void CellBoundaryConnection<IndexType>::setNextBElemWRTCR(IndexType nextIndex){
    nextBElemWRTCR = nextIndex;
}

template<typename IndexType>
void CellBoundaryConnection<IndexType>::setNextBElemWRTCL(IndexType nextIndex){
    nextBElemWRTCL = nextIndex;
}

template<typename IndexType>
bool CellBoundaryConnection<IndexType>::setNextBElem(IndexType nextBElemIndex, IndexType cellIndex){

    // CellIndex is invalid then false returned
    if (cellIndex == INVALID_INDEX(IndexType)){
        return false;
    }
    // first test wether cell index eaqules left or right
    // then is posible to set up invalid indexes
    if (CellConnection<IndexType>::getCellLeftIndex() == cellIndex) {

        setNextBElemWRTCL(nextBElemIndex);
        return true;

    }

    if (CellConnection<IndexType>::getCellRightIndex() == cellIndex){

        setNextBElemWRTCR(nextBElemIndex);
        return true;

    }


    // Attribute CellRightIndex is invalid then
    // CellRightIndex is set as cellIndex
    if (CellConnection<IndexType>::getCellLeftIndex() == INVALID_INDEX(IndexType)){

        CellConnection<IndexType>::setCellLeftIndex(cellIndex);

    }

    // Parameter cellIndex is equal to CellRightIndex
    // then nextEdgeWRTCL is set as nextEdgeIndex, ret true
    if (cellIndex == CellConnection<IndexType>::getCellLeftIndex()) {

        setNextBElemWRTCL(nextBElemIndex);

        return true;

    } else {

        // Attribute CellRightIndex is invalid
        // but CellLeftIndex is already filled
        // then set CellLeftIndex as cellIndex
        if(CellConnection<IndexType>::getCellRightIndex() == INVALID_INDEX(IndexType)){
            CellConnection<IndexType>::setCellRightIndex(cellIndex);
        }

        // Parameter cellIndex is equal to CellLeftIndex then
        // set NextEdgeWRTCR as nextEdgeIndex, ret true
        if (cellIndex == CellConnection<IndexType>::getCellRightIndex()){

            setNextBElemWRTCR(nextBElemIndex);

            return true;



            // If both CellLeftIndex and CellRightIndex
            // has a non -1 value
            // and parameter cellIndex does not
            // equal CellRightIndex neither CellLeftIndex
            // then return false
        } else {

            return false;

        }
    }
}

template<typename IndexType>
IndexType CellBoundaryConnection<IndexType>::getNextBElemWRTCL() const {
    return nextBElemWRTCL;
}

template<typename IndexType>
IndexType CellBoundaryConnection<IndexType>::getNextBElemWRTCR() const {
    return nextBElemWRTCR;
}

template<typename IndexType>
IndexType CellBoundaryConnection<IndexType>::getNextBElem(IndexType cellIndex) const{

    // If cell is invalied then
    if (cellIndex == INVALID_INDEX(IndexType)) {
        throw std::runtime_error("Invalid index given to the getNextBElem");
    }

    // If the cell is equal the Cell1 then return the NextBElemWRTCR
    if(cellIndex == this->getCellRightIndex()){
        return getNextBElemWRTCR();

    // If the cell is equal the Cell2 then return the NextBElemWRTCL
    } else if (cellIndex == this->getCellLeftIndex()){
        return getNextBElemWRTCL();

    // If the cell is not equal left cell neither cell right then return invalid index
    } else {
        std::stringstream error;
        error << "Neither of cell indexes (" << this->getCellLeftIndex() << ","
              << this->getCellRightIndex() << ") matches the given one ("
              << cellIndex << ")";

        throw std::runtime_error(error.str());
    }
}




#endif // CELLBOUNDARYCONNECTION_H
