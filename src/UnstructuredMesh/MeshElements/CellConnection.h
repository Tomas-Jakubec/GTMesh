#ifndef CELLCONNECTION_H
#define CELLCONNECTION_H
#include "../UnstructuredMeshDefine.h"
#include <sstream>
#include <stdexcept>

/**
 * @brief The CellConnection class
*/
template<typename IndexType>
class CellConnection {
    // Indexes to two cells which neighbours
    // with this edge
    IndexType cellRightIndex, cellLeftIndex;
    public:
    /**
     * @brief CellConnection
     * @param cellRight
     * @param cellLeft
     */
    CellConnection( IndexType cellLeft = INVALID_INDEX(IndexType),
                    IndexType cellRight = INVALID_INDEX(IndexType) );

    /**
     * @brief setCellRightIndex
     * @param cellIndex
     */
    void setCellRightIndex( IndexType cellIndex );

    /**
     * @brief SetCellLeftIndex
     * @param cellIndex
     */
    void setCellLeftIndex( IndexType cellIndex );

    /**
     * @brief SetCellIndex
     * @param cellIndex
     * @return
     */
    bool setCellIndex( IndexType cellIndex );

    /**
     * @brief GetCellRightIndex
     * @return
     */
    IndexType getCellRightIndex() const;

    /**
     * @brief GetCellLeftIndex
     * @return
     */
    IndexType getCellLeftIndex() const;

    // Returns the other Cell than sent by parameter
    /**
     * @brief GetOtherCellIndex
     * @param cellIndex
     * @return
     */
    IndexType getOtherCellIndex( IndexType cellIndex ) const;


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

template<typename IndexType>
CellConnection<IndexType>::CellConnection(IndexType cellLeft, IndexType cellRight){
    cellRightIndex = cellRight;
    cellLeftIndex = cellLeft;
}

template<typename IndexType>
void CellConnection<IndexType>::setCellRightIndex(IndexType cellIndex){
    cellRightIndex = cellIndex;
}

template<typename IndexType>
void CellConnection<IndexType>::setCellLeftIndex(IndexType cellIndex){
    cellLeftIndex = cellIndex;
}

template<typename IndexType>
bool CellConnection<IndexType>::setCellIndex(IndexType cellIndex){
    // If the parameter cell is nullptr then ret false
    if (cellIndex == INVALID_INDEX(IndexType)){
        return false;
    }

    // If the CellLeftIndex is lower than 0
    // then set CellLeftIndex as cellIndex, ret true
    if (getCellLeftIndex() == INVALID_INDEX(IndexType)) {

        setCellLeftIndex(cellIndex);

        return true;
        // If the CellRightIndex is valid
        // and CellLeftIndex is not
        // then set CellRightIndex as cellIndex
    } else if (getCellRightIndex() == INVALID_INDEX(IndexType)) {

        setCellRightIndex(cellIndex);

        return true;
        // If both CellLeftIndex and CellRightIndex are >= 0
        // then return false
    } else {

        return false;

    }


}

template<typename IndexType>
IndexType CellConnection<IndexType>::getCellRightIndex() const {
    return cellRightIndex;
}

template<typename IndexType>
IndexType CellConnection<IndexType>::getCellLeftIndex() const {
    return cellLeftIndex;
}

template<typename IndexType>
IndexType CellConnection<IndexType>::getOtherCellIndex(IndexType cellIndex) const{
    // If cell is invalied then
    if (cellIndex == INVALID_INDEX(IndexType)) {
        throw std::runtime_error("Invalid index given to the getOtherCellIndex");
    }

    // If the cell is equal the Cell1 then return the NextBElemWRTCR
    if(cellIndex == this->getCellRightIndex()){
        return getCellLeftIndex();

    // If the cell is equal the Cell2 then return the NextBElemWRTCL
    } else if (cellIndex == this->getCellLeftIndex()){
        return getCellRightIndex();

    // If the cell is not equal left cell neither cell right then return invalid index
    } else {
        std::stringstream error;
        error << "Neither of cell indexes (" << this->getCellLeftIndex() << ","
              << this->getCellRightIndex() << ") matches the given one ("
              << cellIndex << ")";

        throw std::runtime_error(error.str());
    }
}

template<typename IndexType>
bool CellConnection<IndexType>::cellsOK() const {
    return getCellRightIndex() != INVALID_INDEX(IndexType) && getCellLeftIndex() != INVALID_INDEX(IndexType);
}

template<typename IndexType>
void CellConnection<IndexType>::swapCellsLR() {


    IndexType cellLeftIndex = getCellRightIndex();
    setCellRightIndex(getCellLeftIndex());
    setCellLeftIndex(cellLeftIndex);


}
#endif // CELLCONNECTION_H
