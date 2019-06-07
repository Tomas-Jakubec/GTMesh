#ifndef CELLBOUNDARYCONNECTION_H
#define CELLBOUNDARYCONNECTION_H
#include "unstructed_mesh_define.h"
#include "cellconnection.h"

template<typename indexType>
class CellBoundaryConnection : public CellConnection<indexType> {
    /**
     * @brief NextBElemWRTCL
     *
     * Index of the next boundary element with respect to the left cell
     */
    indexType NextBElemWRTCL;
    /**
     * @brief NextBElemWRTCR
     *
     * Index of the next boundary element with respect to right cell
     */
    indexType NextBElemWRTCR;
public:
    CellBoundaryConnection(indexType cellLeft = INVALID_INDEX(indexType),
                           indexType cellRight = INVALID_INDEX(indexType),
                           indexType nextLeft = INVALID_INDEX(indexType),
                           indexType nextRight = INVALID_INDEX(indexType));
    /*
    ** Set atribute methods
    */
    void SetNextBElemWRTCR(indexType nextIndex);

    void SetNextBElemWRTCL(indexType nextIndex);
    bool SetNextEdge(indexType nextEdgeIndex, indexType cellIndex);

    indexType GetNextBElemWRTCL() const;

    indexType GetNextBElemWRTCR() const {
        return NextBElemWRTCR;
    }

    indexType GetNextBElem(indexType cellIndex) const;

    // Returns the other Cell than sent by parameter
    indexType GetOtherCellIndex(indexType cellIndex) const;


};

/**
 * Implementation of methods
 */

template<typename indexType>
CellBoundaryConnection<indexType>::CellBoundaryConnection(indexType cellLeft, indexType cellRight, indexType nextLeft, indexType nextRight)
    : CellConnection<indexType> (cellLeft, cellRight) {
    NextBElemWRTCL = nextLeft;
    NextBElemWRTCR = nextRight;
}

template<typename indexType>
void CellBoundaryConnection<indexType>::SetNextBElemWRTCR(indexType nextIndex){
    NextBElemWRTCR = nextIndex;
}

template<typename indexType>
void CellBoundaryConnection<indexType>::SetNextBElemWRTCL(indexType nextIndex){
    NextBElemWRTCL = nextIndex;
}

template<typename indexType>
bool CellBoundaryConnection<indexType>::SetNextEdge(indexType nextEdgeIndex, indexType cellIndex){

    // CellIndex is invalid then false returned
    if (cellIndex == INVALID_INDEX(indexType)){
        return false;
    }
    // first test wether cell index eaqules left or right
    // then is posible to set up invalid indexes
    if (CellConnection<indexType>::GetCellLeftIndex() == cellIndex) {

        SetNextBElemWRTCL(nextEdgeIndex);
        return true;

    }

    if (CellConnection<indexType>::GetCellRightIndex() == cellIndex){

        SetNextBElemWRTCR(nextEdgeIndex);
        return true;

    }


    // Attribute CellRightIndex is invalid then
    // CellRightIndex is set as cellIndex
    if (CellConnection<indexType>::GetCellLeftIndex() == INVALID_INDEX(indexType)){

        SetCellLeftIndex(cellIndex);

    }

    // Parameter cellIndex is equal to CellRightIndex
    // then nextEdgeWRTCL is set as nextEdgeIndex, ret true
    if (cellIndex == CellConnection<indexType>::GetCellLeftIndex()) {

        SetNextBElemWRTCL(nextEdgeIndex);

        return true;

    } else {

        // Attribute CellRightIndex is invalid
        // but CellLeftIndex is already filled
        // then set CellLeftIndex as cellIndex
        if(CellConnection<indexType>::GetCellRightIndex() == INVALID_INDEX(indexType)){
            SetCellRightIndex(cellIndex);
        }

        // Parameter cellIndex is equal to CellLeftIndex then
        // set NextEdgeWRTCR as nextEdgeIndex, ret true
        if (cellIndex == CellConnection<indexType>::GetCellRightIndex()){

            SetNextBElemWRTCR(nextEdgeIndex);

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
indexType CellBoundaryConnection<indexType>::GetNextBElemWRTCL() const {
    return NextBElemWRTCL;
}

template<typename indexType>
indexType CellBoundaryConnection<indexType>::GetNextBElem(indexType cellIndex) const{
    // If cell is nullptr then ret nullptr
    if (cellIndex == INVALID_INDEX(indexType)) {
        return INVALID_INDEX(indexType);
    }


    // If the cell is equal the Cell1 then return the NextBElemWRTCR
    if(cellIndex == CellConnection<indexType>::GetCellRightIndex()){
        return GetNextBElemWRTCR();

        // If the cell is equal the Cell2 then return the NextBElemWRTCL
    } else if (cellIndex == CellConnection<indexType>::GetCellLeftIndex()){
        return GetNextBElemWRTCL();

        // If the cell is not equal left cell neither cell right then return invalid index
    } else {
        return INVALID_INDEX(indexType);
    }
}

template<typename indexType>
indexType CellBoundaryConnection<indexType>::GetOtherCellIndex(indexType cellIndex) const{
    if (cellIndex == CellConnection<indexType>::GetCellLeftIndex()) {
        return CellConnection<indexType>::GetCellRightIndex();
    } else if (cellIndex == CellConnection<indexType>::GetCellRightIndex()){
        return CellConnection<indexType>::GetCellLeftIndex();
    }
    return INVALID_INDEX(indexType);
}



#endif // CELLBOUNDARYCONNECTION_H
