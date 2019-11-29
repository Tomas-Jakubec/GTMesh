#ifndef CELLSDISTANCE_H
#define CELLSDISTANCE_H
#include "../MeshElements/MeshElement.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "../../Debug/Debug.h"
#include "MeshApply.h"
#include <valarray>
#include <set>
#include <map>






template <unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MeshDataContainer<Real, Dimension-1> ComputeCellsDistance(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

    MeshDataContainer<Real, Dimension-1> distances(mesh);

    if(mesh.getBoundaryCells().empty()){
        for(auto& face : mesh.getFaces()) {
            if (face.getCellLeftIndex() != INVALID_INDEX(IndexType) &&
                face.getCellRightIndex() != INVALID_INDEX(IndexType)){

                distances.at(face) = (mesh.getCells().at(face.getCellLeftIndex()).getCenter() -
                                      mesh.getCells().at(face.getCellRightIndex()).getCenter()).normEukleid();

            } else if(face.getCellLeftIndex() != INVALID_INDEX(IndexType) &&
                      face.getCellRightIndex() == INVALID_INDEX(IndexType)){

                distances.at(face) = (mesh.getCells().at(face.getCellLeftIndex()).getCenter() -
                                      face.getCenter()).normEukleid();

            } else if(face.getCellLeftIndex() == INVALID_INDEX(IndexType) &&
                      face.getCellRightIndex() != INVALID_INDEX(IndexType)){

                distances.at(face) = (mesh.getCells().at(face.getCellRightIndex()).getCenter() -
                                      face.getCenter()).normEukleid();
            }
        }

    } else {

        for(auto& face : mesh.getFaces()) {
            auto& cellLeft = (face.getCellLeftIndex() & BOUNDARY_INDEX(IndexType)) == BOUNDARY_INDEX(IndexType)?
                      mesh.getBoundaryCells().at(face.getCellLeftIndex()&EXTRACTING_INDEX(IndexType)):
                      mesh.getCells().at(face.getCellLeftIndex());
            auto& cellRight = (face.getCellRightIndex() & BOUNDARY_INDEX(IndexType)) == BOUNDARY_INDEX(IndexType)?
                      mesh.getBoundaryCells().at(face.getCellRightIndex()&EXTRACTING_INDEX(IndexType)):
                      mesh.getCells().at(face.getCellRightIndex());

            distances.at(face) = (cellLeft.getCenter() - cellRight.getCenter()).normEukleid();
        }
    }


    return distances;
}



#endif // CELLSDISTANCE_H
