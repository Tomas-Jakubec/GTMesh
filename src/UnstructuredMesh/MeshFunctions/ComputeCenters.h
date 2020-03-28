#ifndef COMPUTECENTER_H
#define COMPUTECENTER_H

#include "MeshFunctionsDefine.h"
#include "MeshApply.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/GrammSchmidt.h"
#include <array>

namespace Impl {



template <unsigned int CurrentDimension, unsigned int MeshDimension, ComputationMethod Method = ComputationMethod::DEFAULT>
struct _ComputeCenters {
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<MeshDimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh) {

        auto& elemCenters = centers.template getDataByDim<CurrentDimension>();
        auto& subElemCenters = centers.template getDataByDim<CurrentDimension - 1>();


        for (IndexType i = 0; i < mesh.template getElements<CurrentDimension>().size(); i++) {

            Real subElemCnt = 0;
            MeshApply<CurrentDimension, CurrentDimension - 1>::apply(
                i,
                mesh,
                [&elemCenters, &subElemCenters, &subElemCnt](IndexType elementIndex, IndexType subelementIndex){
                    elemCenters.at(elementIndex) +=  subElemCenters.at(subelementIndex);
                    subElemCnt++;
                }
            );

            elemCenters.at(i) /= subElemCnt;
        }

        _ComputeCenters<CurrentDimension + 1, MeshDimension, Method>::compute(centers, mesh);
    }
};






template <unsigned int MeshDimension, ComputationMethod Method>
struct _ComputeCenters<MeshDimension, MeshDimension, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Vertex<MeshDimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
                        const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh) {

        auto& elemCenters = centers.template getDataByDim<MeshDimension>();
        auto& subElemCenters = centers.template getDataByDim<MeshDimension - 1>();


        for (IndexType i = 0; i < mesh.template getElements<MeshDimension>().size(); i++) {

            Real subElemCnt = 0;
            MeshApply<MeshDimension, MeshDimension - 1>::apply(
                        i,
                        mesh,
                        [&elemCenters, &subElemCenters, &subElemCnt](IndexType elementIndex, IndexType subelementIndex){

                elemCenters.at(elementIndex) +=  subElemCenters.at(subelementIndex);
                subElemCnt++;
            }
            );

            elemCenters.at(i) /= subElemCnt;
        }
    }

};

template <unsigned int MeshDimension, ComputationMethod Method>
struct _ComputeCenters<1, MeshDimension, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<MeshDimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
            const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){

        std::vector<Vertex<MeshDimension, Real>>& edgeCenters = centers.template getDataByDim<1>();

        for (IndexType edgeIndex = 0; edgeIndex < mesh.template getElements<1>().size(); edgeIndex++) {
            const auto& edge = mesh.getEdges().at(edgeIndex);

            edgeCenters.at(edgeIndex) = (mesh.template getElements<0>().at(edge.getVertexAIndex()) +
                                mesh.template getElements<0>().at(edge.getVertexBIndex())) * 0.5;
        }

        _ComputeCenters<2, MeshDimension, Method>::compute(centers, mesh);
    }
};




/**
 * @brief The _ComputeCenters<_Tp1, _Tp2, CalculationMethod::MESH_TESSELLATED> struct
 * specialization of computation of centers for mesh which needs tessellation of faces
 */
template <>
struct _ComputeCenters<2, 3, ComputationMethod::TESSELLATED> {
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<3, Real>, make_custom_integer_sequence_t<unsigned int, 1, 3>>& centers,
            const MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& elemCenters = centers.template getDataByDim<2>();
        auto& subElemCenters = centers.template getDataByDim<2 - 1>();

        for (IndexType i = 0; i < mesh.template getElements<2>().size(); i++) {

            Real subElemCnt = 0;
            MeshApply<2, 1>::apply(i, mesh,[&](IndexType faceIndex, IndexType edgeIndex){

                elemCenters.at(faceIndex) +=  subElemCenters.at(edgeIndex);
                subElemCnt++;
            });


            elemCenters.at(i) /= subElemCnt;

            Vertex<3, Real> tempVert = {};
            Real surfTotal = 0.0;
            MeshApply<2, 1>::apply(i, mesh,[&](IndexType faceIndex, IndexType edgeIndex){

                IndexType AI = mesh.getEdges().at(edgeIndex).getVertexAIndex();
                IndexType BI = mesh.getEdges().at(edgeIndex).getVertexBIndex();
                std::array<Vertex<3, Real>, 2> v = {elemCenters.at(faceIndex) - mesh.getVertices().at(AI), elemCenters.at(faceIndex) - mesh.getVertices().at(BI)};
                std::array<Real, 2> norms;
                grammSchmidt<2, 3, IndexType, Real>(v, norms);
                Real surf = norms.at(0) * 0.5 * norms.at(1);

                tempVert += subElemCenters.at(edgeIndex) * (surf * (2.0 / 3.0));
                surfTotal += surf;
            });


            elemCenters.at(i) = (elemCenters.at(i) / 3.0) + (tempVert / surfTotal);
        }

        _ComputeCenters<2 + 1, 3, ComputationMethod::TESSELLATED>::compute(centers, mesh);
    }
};

} // namespace Impl

/**
 * @brief computeCenters calculates centers of all elements
 * of higher dimension than 0
 * @param mesh
 */
template <ComputationMethod Method, unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>
computeCenters(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

    MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> centers(mesh);

    Impl::_ComputeCenters<1, Dimension, Method>::compute(centers, mesh);

    return centers;
}


#endif // COMPUTECENTER_H
