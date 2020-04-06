#ifndef COMPUTEMEASURES_H
#define COMPUTEMEASURES_H

#include "MeshFunctionsDefine.h"
#include "MeshApply.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/GramSchmidt.h"
#include <array>



namespace Impl {


template <unsigned int CurrentDimension, unsigned int MeshDimension, ComputationMethod Method = METHOD_DEFAULT>
struct _ComputeMeasures{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>&,MeshElements<MeshDimension, IndexType, Real, Reserve...>&){
        static_assert (MeshDimension <= 3,"The measure computation of mesh of dimension higher than 3 is not implemented yet.");
    }
};


template <unsigned int MeshDimension, ComputationMethod Method>
struct _ComputeMeasures<1, MeshDimension, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& measures,MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){

        auto& edgeLengths = measures.template getDataByDim<1>();

        for (IndexType edgeIndex = 0; edgeIndex < mesh.template getElements<1>().size(); edgeIndex++) {
            auto& edge = mesh.getEdges().at(edgeIndex);
            edgeLengths.at(edgeIndex) = (mesh.getVertices().at(edge.getVertexAIndex()) -
                                         mesh.getVertices().at(edge.getVertexBIndex())).normEuclid();
        }

        _ComputeMeasures<2, MeshDimension>::compute(measures, mesh);
    }
};


template <ComputationMethod Method>
struct _ComputeMeasures<3, 3, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& cellMeasures = measures.template getDataByDim<3>();

        for (IndexType cellIndex = 0; cellIndex < mesh.getCells().size(); cellIndex++) {

            const auto& cell = mesh.getCells().at(cellIndex);
            IndexType tmpFace = cell.getBoundaryElementIndex();
            Real measure = Real();
            const Vertex<3,Real>& cellCenter = cell.getCenter();

            do {
                // select 3 different vertices
                IndexType vAIndex = mesh.getEdges().at(mesh.getFaces().at(tmpFace).getSubelements()[0]).getVertexAIndex();
                IndexType vBIndex = mesh.getEdges().at(mesh.getFaces().at(tmpFace).getSubelements()[0]).getVertexBIndex();

                Vertex<3,Real>& a = mesh.getVertices().at(vAIndex);
                Vertex<3,Real>& b = mesh.getVertices().at(vBIndex);
                Vertex<3,Real>& c = mesh.getFaces().at(tmpFace).getCenter();

                std::array<Vertex<3,Real>, 3> gsVecs = {a - b, a - c, a - cellCenter};
                std::array<Real, 3> gsNorms = {};

                gramSchmidt<3,3,IndexType, Real>(gsVecs, gsNorms);

                Real distance = gsNorms[2];

                Real tmp = distance * measures.template getDataByDim<2>().at(tmpFace);
                measure += tmp / 3.0;

                tmpFace = mesh.getFaces().at(tmpFace).getNextBElem(cellIndex);
            } while (tmpFace != cell.getBoundaryElementIndex());

            cellMeasures.at(cellIndex) = measure;
        }
    }
};

template <ComputationMethod Method>
struct _ComputeMeasures<2, 2, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 2>>& measures,MeshElements<2, IndexType, Real, Reserve...>& mesh){

        auto& surfaceMeasures = measures.template getDataByDim<2>();

        for (IndexType cellIndex = 0; cellIndex < mesh.getCells().size(); cellIndex++) {
            typename MeshElements<2, IndexType, Real, Reserve...>::template ElementType<2>& cell = mesh.getCells().at(cellIndex);
            IndexType tmpEdge = cell.getBoundaryElementIndex();
            Real measure = Real();
            Vertex<2,Real>& cellCenter = cell.getCenter();
            do {
                Vertex<2,Real>& a = mesh.getVertices().at(mesh.getEdges().at(tmpEdge).getVertexAIndex());
                Vertex<2,Real>& b = mesh.getVertices().at(mesh.getEdges().at(tmpEdge).getVertexBIndex());
                double tmp = (cellCenter[0] - a[0]) * (b[1] - a[1]);
                tmp -= (cellCenter[1] - a[1]) * (b[0] - a[0]);
                measure += 0.5 * fabs(tmp);

                tmpEdge = mesh.getEdges().at(tmpEdge).getNextBElem(cellIndex);
            } while (tmpEdge != cell.getBoundaryElementIndex());

            surfaceMeasures.at(cellIndex) = measure;
        }
    }
};



template <ComputationMethod Method>
struct _ComputeMeasures<2, 3, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& surfaceMeasures = measures.template getDataByDim<2>();

        for (IndexType faceIndex = 0; faceIndex < mesh.getFaces().size(); faceIndex++) {

            const auto& face = mesh.template getElements<2>().at(faceIndex);
            Real measure = Real();
            const Vertex<3,Real>& faceCenter = face.getCenter();
            for(auto sube : face.getSubelements()){
                const auto& edge = mesh.getEdges().at(sube);
                Vertex<3,Real>& a = mesh.getVertices().at(edge.getVertexAIndex());
                Vertex<3,Real>& b = mesh.getVertices().at(edge.getVertexBIndex());

                std::array<Vertex<3,Real>, 2> gsVecs = {b - a, faceCenter - a};
                std::array<Real, 2> gsNorms = {};

                gramSchmidt<2,3,IndexType, Real>(gsVecs, gsNorms);

                measure += 0.5 * gsNorms[1] * measures.template getDataByDim<1>().at(sube);
            }
            surfaceMeasures.at(faceIndex) = measure;
        }
        _ComputeMeasures<3, 3>::compute(measures, mesh);
    }
};








template <>
struct _ComputeMeasures<3, 3, METHOD_TESSELLATED>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& cellMeasures = measures.template getDataByDim<3>();

        for (IndexType cellIndex = 0; cellIndex < mesh.getCells().size(); cellIndex++) {
            typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<3>& cell = mesh.getCells().at(cellIndex);
            Vertex<3,Real>& cellCenter = cell.getCenter();
            Real measure = Real();
            MeshApply<3,2>::apply(
                        cellIndex,
                        mesh,
                [&cellCenter, &cellMeasures, &mesh, &measure](IndexType , IndexType faceIndex){

                Vertex<3,Real>& faceCenter = mesh.getFaces().at(faceIndex).getCenter();

                MeshApply<2,1>::apply(faceIndex, mesh, [&](IndexType , IndexType edgeIndex){
                    Vertex<3,Real>& vertA = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexAIndex());
                    Vertex<3,Real>& vertB = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexBIndex());

                    std::array<Vertex<3,Real>, 3> pyramidVec = {vertA - faceCenter, vertB - faceCenter, cellCenter - faceCenter};
                    std::array<Real, 3> norms;
                    gramSchmidt<3, 3, IndexType, Real>(pyramidVec, norms);

                    measure += norms.at(0) * norms.at(1) * norms.at(2) * (1.0/6.0);

                });
                }

            );

            cellMeasures.at(cellIndex) = measure;
        }
    }
};

} // namespace Impl





template <ComputationMethod Method, unsigned int MeshDimension,typename IndexType, typename Real, unsigned int ...Reserve>
MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>> computeMeasures(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>> measures(mesh);

    Impl::_ComputeMeasures<1, MeshDimension, Method>::compute(measures, mesh);

    return measures;
}


#endif // COMPUTEMEASURES_H
