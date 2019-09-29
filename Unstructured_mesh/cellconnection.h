#ifndef CELLCONNECTION_H
#define CELLCONNECTION_H
#include "unstructed_mesh_define.h"


/**
 * @brief The CellConnection class
*/
    template<typename indexType>
    class CellConnection {
        // Indexes to two cells which neighbours
        // with this edge
        indexType cellRightIndex, cellLeftIndex;
    public:
    /**
     * @brief CellConnection
     * @param cellRight
     * @param cellLeft
     */
    CellConnection(indexType cellLeft = INVALID_INDEX(indexType),
                   indexType cellRight = INVALID_INDEX(indexType));

    /**
     * @brief setCellRightIndex
     * @param cellIndex
     */
    void setCellRightIndex(indexType cellIndex);

    /**
     * @brief SetCellLeftIndex
     * @param cellIndex
     */
    void setCellLeftIndex(indexType cellIndex);

    /**
     * @brief SetCellIndex
     * @param cellIndex
     * @return
     */
    bool setCellIndex(indexType cellIndex);

    /**
     * @brief GetCellRightIndex
     * @return
     */
    indexType getCellRightIndex() const;

    /**
     * @brief GetCellLeftIndex
     * @return
     */
    indexType getCellLeftIndex() const;

    // Returns the other Cell than sent by parameter
    /**
     * @brief GetOtherCellIndex
     * @param cellIndex
     * @return
     */
    indexType getOtherCellIndex(indexType cellIndex) const;


    /**
     * @brief CellsOK
     * @return true if both cell indexes are set
     */
    bool cellsOK() const;


    /**
     * @brief SwapCellsLR
     *
     * swaps the left and right cell indexes
     */
    void swapCellsLR();

};


/***
 *
 * Implementation of methods
 *
*/

template<typename indexType>
CellConnection<indexType>::CellConnection(indexType cellLeft, indexType cellRight){
    cellRightIndex = cellRight;
    cellLeftIndex = cellLeft;
}

template<typename indexType>
void CellConnection<indexType>::setCellRightIndex(indexType cellIndex){
    cellRightIndex = cellIndex;
}

template<typename indexType>
void CellConnection<indexType>::setCellLeftIndex(indexType cellIndex){
    cellLeftIndex = cellIndex;
}

template<typename indexType>
bool CellConnection<indexType>::setCellIndex(indexType cellIndex){
    // If the parameter cell is nullptr then ret false
    if (cellIndex == INVALID_INDEX(indexType)){
        return false;
    }

    // If the CellLeftIndex is lower than 0
    // then set CellLeftIndex as cellIndex, ret true
    if (getCellLeftIndex() == INVALID_INDEX(indexType)) {

        setCellLeftIndex(cellIndex);

        return true;
        // If the CellRightIndex is valid
        // and CellLeftIndex is not
        // then set CellRightIndex as cellIndex
    } else if (getCellRightIndex() == INVALID_INDEX(indexType)) {

        setCellRightIndex(cellIndex);

        return true;
        // If both CellLeftIndex and CellRightIndex are >= 0
        // then return false
    } else {

        return false;

    }


}

template<typename indexType>
indexType CellConnection<indexType>::getCellRightIndex() const {
    return cellRightIndex;
}

template<typename indexType>
indexType CellConnection<indexType>::getCellLeftIndex() const {
    return cellLeftIndex;
}

template<typename indexType>
indexType CellConnection<indexType>::getOtherCellIndex(indexType cellIndex) const{
    if (cellIndex == getCellLeftIndex()) {
        return getCellRightIndex();
    } else if (cellIndex == getCellRightIndex()){
        return getCellLeftIndex();
    }
    return INVALID_INDEX(indexType);
}

template<typename indexType>
bool CellConnection<indexType>::cellsOK() const {
    return getCellRightIndex() != INVALID_INDEX(indexType) && getCellLeftIndex() != INVALID_INDEX(indexType);
}

template<typename indexType>
void CellConnection<indexType>::swapCellsLR() {


    indexType cellLeftIndex = getCellRightIndex();
    setCellRightIndex(getCellLeftIndex());
    setCellLeftIndex(cellLeftIndex);


}
#endif // CELLCONNECTION_H
