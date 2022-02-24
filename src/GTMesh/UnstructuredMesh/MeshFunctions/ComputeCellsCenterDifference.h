#ifndef COMPUTECENTERDIFFERENCE_H
#define COMPUTECENTERDIFFERENCE_H
#include "../../NumericStaticArray/Vector.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../MeshElements/MeshElements.h"
#include "MeshApply.h"
#include <map>
#include <set>
#include <valarray>

template<unsigned int Dimension, typename IndexType, typename Real, unsigned int... Reserve>
MeshDataContainer<Vector<Dimension, Real>, Dimension - 1> computeCellsCenterDifference(
    const MeshElements<Dimension, IndexType, Real, Reserve...> &mesh)
{
    MeshDataContainer<Vector<Dimension, Real>, Dimension - 1> connectingLines(mesh);

    if (mesh.getBoundaryCells().empty()) {
        for (auto &face : mesh.getFaces()) {
            if (!isInvalidIndex(face.getCellLeftIndex()) && !isInvalidIndex(face.getCellRightIndex())) {
                connectingLines.at(face) = mesh.getCells().at(face.getCellRightIndex()).getCenter()
                                           - mesh.getCells().at(face.getCellLeftIndex()).getCenter();

            } else if (!isInvalidIndex(face.getCellLeftIndex()) && isInvalidIndex(face.getCellRightIndex())) {
                connectingLines.at(face) = face.getCenter() - mesh.getCells().at(face.getCellLeftIndex()).getCenter();

            } else if (isInvalidIndex(face.getCellLeftIndex()) && !isInvalidIndex(face.getCellRightIndex())) {
                connectingLines.at(face) = face.getCenter() - mesh.getCells().at(face.getCellRightIndex()).getCenter();
            }
        }

    } else {
        for (auto &face : mesh.getFaces()) {
            auto &cellLeft = isBoundaryIndex(face.getCellLeftIndex())
                                 ? mesh.getBoundaryCells().at(extractBoundaryIndex(face.getCellLeftIndex()))
                                 : mesh.getCells().at(face.getCellLeftIndex());
            auto &cellRight = isBoundaryIndex(face.getCellRightIndex())
                                  ? mesh.getBoundaryCells().at(extractBoundaryIndex(face.getCellRightIndex()))
                                  : mesh.getCells().at(face.getCellRightIndex());

            connectingLines.at(face) = cellRight.getCenter() - cellLeft.getCenter();
        }
    }

    return connectingLines;
}

#endif // COMPUTECENTERDIFFERENCE_H
