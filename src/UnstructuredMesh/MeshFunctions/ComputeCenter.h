#ifndef COMPUTECENTER_H
#define COMPUTECENTER_H

#include "MeshFunctionsDefine.h"
#include "MeshApply.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "GrammSchmidt.h"




template <unsigned int dim, unsigned int Dimension, ComputationMethod Method = ComputationMethod::DEFAULT>
struct _ComputeCenters {
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& centers,
            MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& elemCenters = centers.template getDataByDim<dim>();
        auto& subElemCenters = centers.template getDataByDim<dim - 1>();


        for (IndexType i = 0; i < mesh.template getElements<dim>().size(); i++) {
            auto& element = mesh.template getElements<dim>().at(i);

            Real subElemCnt = 0;
            for(auto& sub : element.getSubelements()){
                elemCenters.at(i) +=  subElemCenters.at(sub.index);
                subElemCnt++;
            }

            elemCenters.at(i) /= subElemCnt;
        }

        _ComputeCenters<dim + 1, Dimension, Method>::compute(centers, mesh);
    }
};






template <unsigned int Dimension, ComputationMethod Method>
struct _ComputeCenters<Dimension, Dimension, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& elemCenters = centers.template getDataByDim<Dimension>();
        auto& subElemCenters = centers.template getDataByDim<Dimension - 1>();


        for (IndexType i = 0; i < mesh.template getElements<Dimension>().size(); i++) {
            auto& element = mesh.template getElements<Dimension>().at(i);

            Real subElemCnt = 0;
            IndexType tmpFaceIndex = element.getBoundaryElementIndex();
            do {
                elemCenters.at(i) +=  subElemCenters.at(tmpFaceIndex);
                subElemCnt++;
                tmpFaceIndex = mesh.getFaces()[tmpFaceIndex].getNextBElem(i);
            } while (tmpFaceIndex != element.getBoundaryElementIndex());

            elemCenters.at(i) /= subElemCnt;
        }
    }

};

template <unsigned int Dimension, ComputationMethod Method>
struct _ComputeCenters<1, Dimension, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& centers,
            MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        std::vector<Vertex<Dimension, Real>>& edgeCenters = centers.template getDataByDim<1>();

        for (auto& edge : mesh.template getElements<1>()) {

            edgeCenters.at(edge.getIndex()) = (mesh.template getElements<0>().at(edge.getVertexAIndex()) +
                                mesh.template getElements<0>().at(edge.getVertexBIndex())) * 0.5;
        }

        _ComputeCenters<2, Dimension, Method>::compute(centers, mesh);
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
            MeshElements<3, IndexType, Real, Reserve...>& mesh){
        DBGMSG("centers computation tessellated");
        auto& elemCenters = centers.template getDataByDim<2>();
        auto& subElemCenters = centers.template getDataByDim<2 - 1>();

        for (IndexType i = 0; i < mesh.template getElements<2>().size(); i++) {
            auto& element = mesh.template getElements<2>().at(i);

            Real subElemCnt = 0;
            for(auto& sub : element.getSubelements()){
                elemCenters.at(i) +=  subElemCenters.at(sub.index);
                subElemCnt++;
            }

            elemCenters.at(i) /= subElemCnt;
        }

        // calculate volume
        for (IndexType i = 0; i < mesh.template getElements<2>().size(); i++) {
            auto& element = mesh.template getElements<2>().at(i);

            Vertex<3, Real> tempVert = {};
            Real surfTotal = 0.0;
            for(auto& sub : element.getSubelements()){
                IndexType AI = mesh.getEdges().at(sub.index).getVertexAIndex();
                IndexType BI = mesh.getEdges().at(sub.index).getVertexBIndex();
                std::vector<Vertex<3, Real>> v = {elemCenters.at(i) - mesh.getVertices().at(AI), elemCenters.at(i) - mesh.getVertices().at(BI)};

                GrammSchmidt<3, IndexType, Real>(v);
                Real surf = v.at(0).normEukleid() * 0.5 * v.at(1).normEukleid();

                tempVert += subElemCenters.at(sub.index) * (surf * (2.0 / 3.0));
                surfTotal += surf;
            }

            elemCenters.at(i) = (elemCenters.at(i) / 3.0) + (tempVert / surfTotal);
        }

        _ComputeCenters<2 + 1, 3, ComputationMethod::TESSELLATED>::compute(centers, mesh);
    }
};



template <ComputationMethod Method, unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>
ComputeCenters(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

     MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> centers(mesh);

    _ComputeCenters<1, Dimension, Method>::compute(centers, mesh);

    return centers;
}


#endif // COMPUTECENTER_H
