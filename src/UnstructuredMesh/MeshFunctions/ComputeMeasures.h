#ifndef COMPUTEMEASURES_H
#define COMPUTEMEASURES_H

#include "MeshFunctionsDefine.h"
#include "MeshApply.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/GrammSchmidt.h"
#include <array>



namespace Detail {


template <unsigned int dim, unsigned int Dimension, ComputationMethod Method = DEFAULT>
struct _ComputeMeasures{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>&,MeshElements<Dimension, IndexType, Real, Reserve...>&){
        static_assert (Dimension <= 3,"The measure computation of mesh of dimension higher than 3 is not implemented yet.");
        throw std::runtime_error("The measure computation of mesh of dimension higher than 3 is not implemented yet.");
    }
};


template <unsigned int Dimension, ComputationMethod Method>
struct _ComputeMeasures<1, Dimension, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& measures,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& edgeLengths = measures.template getDataByDim<1>();

        for (auto& edge : mesh.getEdges()) {
            edgeLengths.at(edge.getIndex()) = (mesh.getVertices().at(edge.getVertexAIndex()) -
                                               mesh.getVertices().at(edge.getVertexBIndex())).normEukleid();
        }

        _ComputeMeasures<2, Dimension>::compute(measures, mesh);
    }
};


template <ComputationMethod Method>
struct _ComputeMeasures<3, 3, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& cellMeasures = measures.template getDataByDim<3>();

        for (typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<3>& cell : mesh.getCells()) {
            IndexType tmpFace = cell.getBoundaryElementIndex();
            Real measure = Real();
            Vertex<3,Real>& cellCenter = cell.getCenter();

            do {
                // select 3 different vertices
                IndexType vAIndex = mesh.getEdges().at(mesh.getFaces().at(tmpFace).getSubelements()[0].index).getVertexAIndex();
                IndexType vBIndex = mesh.getEdges().at(mesh.getFaces().at(tmpFace).getSubelements()[0].index).getVertexBIndex();
                IndexType vCIndex = mesh.getEdges().at(mesh.getFaces().at(tmpFace).getSubelements()[1].index).getVertexAIndex();
                if(vCIndex == vAIndex || vCIndex == vBIndex) {
                    vCIndex = mesh.getEdges().at(mesh.getFaces().at(tmpFace).getSubelements()[1].index).getVertexBIndex();
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

                Real distance = (vAmcC + (vBmA * param_t) + (vCmA * param_s)).normEukleid();

                Real tmp = distance * measures.template getDataByDim<2>().at(tmpFace);
                measure += tmp / 3.0;

                tmpFace = mesh.getFaces().at(tmpFace).getNextBElem(cell.getIndex());
            } while (tmpFace != cell.getBoundaryElementIndex());

            cellMeasures.at(cell.getIndex()) = measure;
        }
    }
};

template <ComputationMethod Method>
struct _ComputeMeasures<2, 2, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 2>>& measures,MeshElements<2, IndexType, Real, Reserve...>& mesh){

        auto& surfaceMeasures = measures.template getDataByDim<2>();

        for (typename MeshElements<2, IndexType, Real, Reserve...>::template ElementType<2>& cell : mesh.getCells()) {
            IndexType tmpEdge = cell.getBoundaryElementIndex();
            Real measure = Real();
            Vertex<2,Real>& cellCenter = cell.getCenter();
            do {
                Vertex<2,Real>& a = mesh.getVertices().at(mesh.getEdges().at(tmpEdge).getVertexAIndex());
                Vertex<2,Real>& b = mesh.getVertices().at(mesh.getEdges().at(tmpEdge).getVertexBIndex());
                double tmp = (cellCenter[0] - a[0]) * (b[1] - a[1]);
                tmp -= (cellCenter[1] - a[1]) * (b[0] - a[0]);
                measure += 0.5 * fabs(tmp);

                tmpEdge = mesh.getEdges().at(tmpEdge).getNextBElem(cell.getIndex());
            } while (tmpEdge != cell.getBoundaryElementIndex());

            surfaceMeasures.at(cell.getIndex()) = measure;
        }
    }
};



template <ComputationMethod Method>
struct _ComputeMeasures<2, 3, Method>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& surfaceMeasures = measures.template getDataByDim<2>();

        for (typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<2>& face : mesh.template getElements<2>()) {

            Real measure = Real();
            Vertex<3,Real>& faceCenter = face.getCenter();
            for(auto sube : face.getSubelements()){

                Vertex<3,Real>& a = mesh.getVertices().at(mesh.getEdges().at(sube.index).getVertexAIndex());
                Vertex<3,Real>& b = mesh.getVertices().at(mesh.getEdges().at(sube.index).getVertexBIndex());

                Real distance = Real();

                Real param = -1.0*(((a-faceCenter)*(b-a))/((b-a).sumOfSquares()));

                distance = (a-faceCenter+(b-a)*param).normEukleid();

                Real tmp = distance * measures.template getDataByDim<1>().at(sube.index);
                measure += tmp * 0.5;
            }
            surfaceMeasures.at(face.getIndex()) = measure;
        }
        _ComputeMeasures<3, 3>::compute(measures, mesh);
    }
};








template <>
struct _ComputeMeasures<3, 3, TESSELLATED>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& cellMeasures = measures.template getDataByDim<3>();

        for (typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<3>& cell : mesh.getCells()) {
            Vertex<3,Real>& cellCenter = cell.getCenter();
            Real measure = Real();
            MeshApply<3,2,3>::apply(
                        cell.getIndex(),
                        mesh,
                [&cellCenter, &cellMeasures, &mesh, &measure](IndexType , IndexType faceIndex){

                Vertex<3,Real>& faceCenter = mesh.getFaces().at(faceIndex).getCenter();

                MeshApply<2,1,3>::apply(faceIndex, mesh, [&](IndexType , IndexType edgeIndex){
                    Vertex<3,Real>& vertA = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexAIndex());
                    Vertex<3,Real>& vertB = mesh.getVertices().at(mesh.getEdges().at(edgeIndex).getVertexBIndex());

                    std::array<Vertex<3,Real>, 3> pyramidVec = {vertA - faceCenter, vertB - faceCenter, cellCenter - faceCenter};
                    std::array<Real, 3> norms;
                    grammSchmidt<3, 3, IndexType, Real>(pyramidVec, norms);

                    measure += norms.at(0) * norms.at(1) * norms.at(2) * (1.0/6.0);

                });
                }

            );

            cellMeasures.at(cell.getIndex()) = measure;
        }
    }
};

} // namespace Detail





template <ComputationMethod Method, unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> ComputeMeasures(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> measures(mesh);

    Detail::_ComputeMeasures<1, Dimension, Method>::compute(measures, mesh);

    return measures;
}


#endif // COMPUTEMEASURES_H
