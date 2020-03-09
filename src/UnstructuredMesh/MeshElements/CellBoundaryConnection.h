#ifndef CELLBOUNDARYCONNECTION_H
#define CELLBOUNDARYCONNECTION_H
#include "../UnstructuredMeshDefine.h"
#include "CellConnection.h"

template<typename indexType>
class CellBoundaryConnection : public CellConnection<indexType> {
    /**
     * @brief nextBElemWRTCL
     *
     * Index of the next boundary element with respect to the left cell
     */
    indexType nextBElemWRTCL;
    /**
     * @brief nextBElemWRTCR
     *
     * Index of the next boundary element with respect to right cell
     */
    indexType nextBElemWRTCR;
public:
    CellBoundaryConnection(indexType cellLeft = INVALID_INDEX(indexType),
                           indexType cellRight = INVALID_INDEX(indexType),
                           indexType nextLeft = INVALID_INDEX(indexType),
                           indexType nextRight = INVALID_INDEX(indexType));
    /*
    ** Set atribute methods
    */
    void setNextBElemWRTCR(indexType nextIndex);

    void setNextBElemWRTCL(indexType nextIndex);

    bool setNextBElem(indexType nextBElemIndex, indexType cellIndex);

    indexType getNextBElemWRTCL() const;

    indexType getNextBElemWRTCR() const;

    indexType getNextBElem(indexType cellIndex) const;

    // Returns the other Cell than sent by parameter
    indexType getOtherCellIndex(indexType cellIndex) const;


};

/**
 * Implementation of methods
 */

template<typename indexType>
CellBoundaryConnection<indexType>::CellBoundaryConnection(indexType cellLeft, indexType cellRight, indexType nextLeft, indexType nextRight)
    : CellConnection<indexType> (cellLeft, cellRight) {
    nextBElemWRTCL = nextLeft;
    nextBElemWRTCR = nextRight;
}

template<typename indexType>
void CellBoundaryConnection<indexType>::setNextBElemWRTCR(indexType nextIndex){
    nextBElemWRTCR = nextIndex;
}

template<typename indexType>
void CellBoundaryConnection<indexType>::setNextBElemWRTCL(indexType nextIndex){
    nextBElemWRTCL = nextIndex;
}

template<typename indexType>
bool CellBoundaryConnection<indexType>::setNextBElem(indexType nextBElemIndex, indexType cellIndex){

    // CellIndex is invalid then false returned
    if (cellIndex == INVALID_INDEX(indexType)){
        return false;
    }
    // first test wether cell index eaqules left or right
    // then is posible to set up invalid indexes
    if (CellConnection<indexType>::getCellLeftIndex() == cellIndex) {

        setNextBElemWRTCL(nextBElemIndex);
        return true;

    }

    if (CellConnection<indexType>::getCellRightIndex() == cellIndex){

        setNextBElemWRTCR(nextBElemIndex);
        return true;

    }


    // Attribute CellRightIndex is invalid then
    // CellRightIndex is set as cellIndex
    if (CellConnection<indexType>::getCellLeftIndex() == INVALID_INDEX(indexType)){

        CellConnection<indexType>::setCellLeftIndex(cellIndex);

    }

    // Parameter cellIndex is equal to CellRightIndex
    // then nextEdgeWRTCL is set as nextEdgeIndex, ret true
    if (cellIndex == CellConnection<indexType>::getCellLeftIndex()) {

        setNextBElemWRTCL(nextBElemIndex);

        return true;

    } else {

        // Attribute CellRightIndex is invalid
        // but CellLeftIndex is already filled
        // then set CellLeftIndex as cellIndex
        if(CellConnection<indexType>::getCellRightIndex() == INVALID_INDEX(indexType)){
            CellConnection<indexType>::setCellRightIndex(cellIndex);
        }

        // Parameter cellIndex is equal to CellLeftIndex then
        // set NextEdgeWRTCR as nextEdgeIndex, ret true
        if (cellIndex == CellConnection<indexType>::getCellRightIndex()){

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

template<typename indexType>
indexType CellBoundaryConnection<indexType>::getNextBElemWRTCL() const {
    return nextBElemWRTCL;
}

template<typename indexType>
indexType CellBoundaryConnection<indexType>::getNextBElemWRTCR() const {
    return nextBElemWRTCR;
}

template<typename indexType>
indexType CellBoundaryConnection<indexType>::getNextBElem(indexType cellIndex) const{
    // If cell is nullptr then ret nullptr
    if (cellIndex == INVALID_INDEX(indexType)) {
        return INVALID_INDEX(indexType);
    }


    // If the cell is equal the Cell1 then return the NextBElemWRTCR
    if(cellIndex == CellConnection<indexType>::getCellRightIndex()){
        return getNextBElemWRTCR();

        // If the cell is equal the Cell2 then return the NextBElemWRTCL
    } else if (cellIndex == CellConnection<indexType>::getCellLeftIndex()){
        return getNextBElemWRTCL();

        // If the cell is not equal left cell neither cell right then return invalid index
    } else {
        return INVALID_INDEX(indexType);
    }
}

template<typename indexType>
indexType CellBoundaryConnection<indexType>::getOtherCellIndex(indexType cellIndex) const{
    if (cellIndex == CellConnection<indexType>::getCellLeftIndex()) {
        return CellConnection<indexType>::getCellRightIndex();
    } else if (cellIndex == CellConnection<indexType>::getCellRightIndex()){
        return CellConnection<indexType>::getCellLeftIndex();
    }
    return INVALID_INDEX(indexType);
}



#endif // CELLBOUNDARYCONNECTION_H
