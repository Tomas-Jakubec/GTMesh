#ifndef COMPUTENORMALS_H
#define COMPUTENORMALS_H

#include "../MeshElements/MeshElement.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "MeshApply.h"


template <unsigned int Dimension>
struct _ComputeNormals{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<Dimension, Real>, Dimension-1>&,MeshElements<Dimension, IndexType, Real, Reserve...>&){
        static_assert (Dimension > 3,"The measure computation of mesh of dimension higher than 3 is not implemented yet.");
        throw std::runtime_error("The computation of face normal vectors of mesh of dimension higher than 3 is not implemented yet.");
    }
};





template <>
struct _ComputeNormals<2>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<2, Real>, 1>& normals,MeshElements<2, IndexType, Real, Reserve...>& mesh){
        for (auto& face : mesh.getEdges()) {
            Vertex<2,Real> a = mesh.getVertices().at(face.getVertexAIndex());
            Vertex<2,Real> b = mesh.getVertices().at(face.getVertexBIndex());
            Vertex<2,Real> dif = b-a;
            normals[face][0] = dif[1];
            normals[face][1] = -dif[0];
            normals[face] /= dif.normEukleid();
        }
    }
};





template <>
struct _ComputeNormals<3>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<3, Real>, 2>& normals,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        for (auto& face : mesh.getFaces()) {

            bool vectorSign = true;
            IndexType cellIndex = face.getCellLeftIndex();
            if (
                    cellIndex == INVALID_INDEX(IndexType) ||
                    (cellIndex & BOUNDARY_INDEX(IndexType)) == BOUNDARY_INDEX(IndexType)
                ) {
                vectorSign = false;
                cellIndex = face.getCellRightIndex();
            }

            Vertex<3,Real>& cellCenter = mesh.getCells().at(cellIndex).getCenter();


            // select 3 different vertices
            IndexType vAIndex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexAIndex();
            IndexType vBIndex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexBIndex();
            IndexType vCIndex = mesh.getEdges().at(face.getSubelements()[1].index).getVertexAIndex();
            if(vCIndex == vAIndex || vCIndex == vBIndex) {
                vCIndex = mesh.getEdges().at(face.getSubelements()[1].index).getVertexBIndex();
            }

            Vertex<3,Real>& a = mesh.getVertices().at(vAIndex);
            Vertex<3,Real>& b = mesh.getVertices().at(vBIndex);
            Vertex<3,Real>& c = mesh.getVertices().at(vCIndex);

            // preparing quiantities
            Vertex<3,Real> vAmcC = (a-cellCenter);
            Vertex<3,Real> vBmA = (b-a);
            Vertex<3,Real> vCmA = (c-a);
            Real inv_sqrBmA = 1.0 / vBmA.sumOfSquares();
            Real inv_sqrCmA = 1.0 / vCmA.sumOfSquares();

            Real denominator = 1.0 / (1.0 - (pow(vCmA*vBmA,2) * inv_sqrBmA * inv_sqrCmA));


            Real param_t = -denominator * (((vAmcC*vBmA) * inv_sqrBmA) - (inv_sqrBmA*inv_sqrCmA*(vAmcC * vCmA)*(vCmA*vBmA)));
            //param_t *= inv_sqrBmA;
            Real param_s = -denominator * (((vAmcC*vCmA) * inv_sqrCmA) - (inv_sqrBmA*inv_sqrCmA*(vAmcC * vBmA)*(vCmA*vBmA)));

            Vertex<3, Real> faceNormal = vAmcC + (vBmA * param_t) + (vCmA * param_s);
            faceNormal /= faceNormal.normEukleid();

            if (!vectorSign) {
                faceNormal *= -1;
            }

            normals.at(face)[0] = fabs(faceNormal[0]) < 1e-8 ? 0 : faceNormal[0];
            normals.at(face)[1] = fabs(faceNormal[1]) < 1e-8 ? 0 : faceNormal[1];
            normals.at(face)[2] = fabs(faceNormal[2]) < 1e-8 ? 0 : faceNormal[2];
        }
    }
};






template <unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MeshDataContainer<Vector<Dimension, Real>, Dimension-1> ComputeFaceNormals(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

    MeshDataContainer<Vector<Dimension, Real>, Dimension-1> normals(mesh);

    _ComputeNormals<Dimension>::compute(normals, mesh);

    return normals;
}


#endif // COMPUTENORMALS_H
