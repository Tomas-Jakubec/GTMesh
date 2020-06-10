#ifndef COMPUTENORMALS_H
#define COMPUTENORMALS_H

#include "../MeshElements/MeshElements.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "../../NumericStaticArray/GramSchmidt.h"
#include "MeshApply.h"
#include "MeshFunctionsDefine.h"

namespace Impl {

template <unsigned int Dimension, ComputationMethod Method = METHOD_DEFAULT>
struct _ComputeNormals{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<Dimension, Real>, Dimension-1>&,
                        const MeshElements<Dimension, IndexType, Real, Reserve...>&){
        static_assert (Dimension > 3,"The measure computation of mesh of dimension higher than 3 is not implemented yet.");
        throw std::runtime_error("The computation of face normal vectors of mesh of dimension higher than 3 is not implemented yet.");
    }
};





template <ComputationMethod Method>
struct _ComputeNormals<2, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<2, Real>, 1>& normals,
                        const MeshElements<2, IndexType, Real, Reserve...>& mesh){
        for (auto& face : mesh.getEdges()) {
            const Vertex<2,Real> a = mesh.getVertices().at(face.getVertexAIndex());
            const Vertex<2,Real> b = mesh.getVertices().at(face.getVertexBIndex());
            Vertex<2,Real> dif = b-a;
            normals[face][0] = dif[1];
            normals[face][1] = -dif[0];
            normals[face] /= dif.normEuclid();
        }
    }
};





template <ComputationMethod Method>
struct _ComputeNormals<3, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<3, Real>, 2>& normals,
                        const MeshElements<3, IndexType, Real, Reserve...>& mesh){
        for (auto& face : mesh.getFaces()) {

            double vectorSign = 1.0;
            IndexType cellIndex = face.getCellLeftIndex();
            if (isInvalidIndex( cellIndex ) || isBoundaryIndex( cellIndex )) {
                vectorSign = -1.0;
                cellIndex = face.getCellRightIndex();
            }

            const Vertex<3,Real>& cellCenter = mesh.getCells().at(cellIndex).getCenter();

            // select 3 different vertices
            IndexType vAIndex = mesh.getEdges().at(face.getSubelements()[0]).getVertexAIndex();
            IndexType vBIndex = mesh.getEdges().at(face.getSubelements()[0]).getVertexBIndex();

            const Vertex<3,Real>& a = mesh.getVertices().at(vAIndex);
            const Vertex<3,Real>& b = mesh.getVertices().at(vBIndex);
            const Vertex<3,Real>& c = face.getCenter();

            std::array<Vertex<3,Real>, 3> gsVecs = {c - a, c - b, c - cellCenter};
            std::array<Real, 3> gsNorms = {};

            gramSchmidt<3,3,IndexType, Real>(gsVecs, gsNorms);

            normals[face] = vectorSign * gsVecs[2];
        }
    }
};


template <>
struct _ComputeNormals<3, METHOD_TESSELLATED>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vector<3, Real>, 2>& normals,
                        const MeshElements<3, IndexType, Real, Reserve...>& mesh){
        for (IndexType faceIndex = 0; faceIndex < mesh.getFaces().size(); faceIndex++) {

            auto& face = mesh.getFaces()[faceIndex];
            bool vectorSign = true;
            IndexType cellIndex = face.getCellLeftIndex();
            if (isInvalidIndex( cellIndex ) || isBoundaryIndex( cellIndex )) {
                vectorSign = false;
                cellIndex = face.getCellRightIndex();
            }

            const Vertex<3, Real>& cellCenter = mesh.getCells().at(cellIndex).getCenter();
            const Vertex<3, Real>& faceCenter = face.getCenter();

            Vector<3, Real> faceNormal = {};

            Real surfTotal = Real();

            MeshApply<2,1>::apply(face.getIndex(), mesh, [&](IndexType , IndexType edgeIndex){
                const Vertex<3,Real>& vertA = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexAIndex());
                const Vertex<3,Real>& vertB = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexBIndex());

                std::array<Vertex<3,Real>, 3> pyramidVec = {faceCenter - vertA, faceCenter - vertB, faceCenter - cellCenter};
                std::array<Real, 3> norms;

                gramSchmidt<3, 3, IndexType, Real>(pyramidVec, norms);

                faceNormal += pyramidVec[2] * 0.5 * norms[0] * norms[1];

                surfTotal += 0.5 * norms[0] * norms[1];

            });

            faceNormal /= surfTotal;

            if (!vectorSign) {
                faceNormal *= -1;
            }
            normals[face] =  faceNormal;
        }
    }
};

}




template <ComputationMethod Method, unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MeshDataContainer<Vector<Dimension, Real>, Dimension-1> computeFaceNormals(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

    MeshDataContainer<Vector<Dimension, Real>, Dimension-1> normals(mesh);

    Impl::_ComputeNormals<Dimension, Method>::compute(normals, mesh);

    return normals;
}


#endif // COMPUTENORMALS_H
