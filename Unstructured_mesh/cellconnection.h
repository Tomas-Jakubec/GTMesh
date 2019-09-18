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
        indexType CellRightIndex, CellLeftIndex;
    public:
    /**
     * @brief CellConnection
     * @param cellRight
     * @param cellLeft
     */
    CellConnection(indexType cellLeft = INVALID_INDEX(indexType),
                   indexType cellRight = INVALID_INDEX(indexType));

    /**
     * @brief SetCellRightIndex
     * @param cellIndex
     */
    void SetCellRightIndex(indexType cellIndex);

    /**
     * @brief SetCellLeftIndex
     * @param cellIndex
     */
    void SetCellLeftIndex(indexType cellIndex);

    /**
     * @brief SetCellIndex
     * @param cellIndex
     * @return
     */
    bool SetCellIndex(indexType cellIndex);

    /**
     * @brief GetCellRightIndex
     * @return
     */
    indexType GetCellRightIndex() const;

    /**
     * @brief GetCellLeftIndex
     * @return
     */
    indexType GetCellLeftIndex() const;

    // Returns the other Cell than sent by parameter
    /**
     * @brief GetOtherCellIndex
     * @param cellIndex
     * @return
     */
    indexType GetOtherCellIndex(indexType cellIndex) const;


    /**
     * @brief CellsOK
     * @return true if both cell indexes are set
     */
    bool CellsOK() const;


    /**
     * @brief SwapCellsLR
     *
     * swaps the left and right cell indexes
     */
    void SwapCellsLR();

};


/***
 *
 * Implementation of methods
 *
*/

template<typename indexType>
CellConnection<indexType>::CellConnection(indexType cellLeft, indexType cellRight){
    CellRightIndex = cellRight;
    CellLeftIndex = cellLeft;
}

template<typename indexType>
void CellConnection<indexType>::SetCellRightIndex(indexType cellIndex){
    CellRightIndex = cellIndex;
}

template<typename indexType>
void CellConnection<indexType>::SetCellLeftIndex(indexType cellIndex){
    CellLeftIndex = cellIndex;
}

template<typename indexType>
bool CellConnection<indexType>::SetCellIndex(indexType cellIndex){
    // If the parameter cell is nullptr then ret false
    if (cellIndex == INVALID_INDEX(indexType)){
        return false;
    }

    // If the CellLeftIndex is lower than 0
    // then set CellLeftIndex as cellIndex, ret true
    if (GetCellLeftIndex() == INVALID_INDEX(indexType)) {

        SetCellLeftIndex(cellIndex);

        return true;
        // If the CellRightIndex is valid
        // and CellLeftIndex is not
        // then set CellRightIndex as cellIndex
    } else if (GetCellRightIndex() == INVALID_INDEX(indexType)) {

        SetCellRightIndex(cellIndex);

        return true;
        // If both CellLeftIndex and CellRightIndex are >= 0
        // then return false
    } else {

        return false;

    }


}

template<typename indexType>
indexType CellConnection<indexType>::GetCellRightIndex() const {
    return CellRightIndex;
}

template<typename indexType>
indexType CellConnection<indexType>::GetCellLeftIndex() const {
    return CellLeftIndex;
}

template<typename indexType>
indexType CellConnection<indexType>::GetOtherCellIndex(indexType cellIndex) const{
    if (cellIndex == GetCellLeftIndex()) {
        return GetCellRightIndex();
    } else if (cellIndex == GetCellRightIndex()){
        return GetCellLeftIndex();
    }
    return INVALID_INDEX(indexType);
}

template<typename indexType>
bool CellConnection<indexType>::CellsOK() const {
    return GetCellRightIndex() != INVALID_INDEX(indexType) && GetCellLeftIndex() != INVALID_INDEX(indexType);
}

template<typename indexType>
void CellConnection<indexType>::SwapCellsLR() {


    indexType cellLeftIndex = GetCellRightIndex();
    SetCellRightIndex(GetCellLeftIndex());
    SetCellLeftIndex(cellLeftIndex);


}
#endif // CELLCONNECTION_H
