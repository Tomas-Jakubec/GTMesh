
#ifndef MESH_FUNCTIONS_H
#define MESH_FUNCTIONS_H
#include "../MeshElements/MeshElement.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "../../Debug/Debug.h"
#include <valarray>
#include <functional>
#include <set>
#include <map>

template <typename Type, Type startIndex, Type EndIndex, int increment = 1, Type... t>
struct MakeCustomIntegerSequence : public MakeCustomIntegerSequence<Type, startIndex + increment, EndIndex, increment, t..., startIndex> {
};

template <typename Type, Type EndIndex, int increment, Type... t>
struct MakeCustomIntegerSequence<Type, EndIndex, EndIndex, increment, t...> {
    using type = std::integer_sequence<Type, t..., EndIndex>;
};

template<typename Type, Type startIndex, Type EndIndex, int increment = 1>
using make_custom_integer_sequence_t = typename MakeCustomIntegerSequence<Type, startIndex, EndIndex, increment>::type;









template <unsigned int dim, unsigned int Dimension>
struct _ComputeCenters{
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

        _ComputeCenters<dim + 1, Dimension>::compute(centers, mesh);
    }
};

template <unsigned int Dimension>
struct _ComputeCenters<Dimension, Dimension>{
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

template <unsigned int Dimension>
struct _ComputeCenters<1, Dimension>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& centers,
            MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        std::vector<Vertex<Dimension, Real>>& edgeCenters = centers.template getDataByDim<1>();

        for (auto& edge : mesh.template getElements<1>()) {

            edgeCenters.at(edge.getIndex()) = (mesh.template getElements<0>().at(edge.getVertexAIndex()) +
                                mesh.template getElements<0>().at(edge.getVertexBIndex())) * 0.5;
        }

        _ComputeCenters<2, Dimension>::compute(centers, mesh);
    }
};




template <unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>
ComputeCenters(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

     MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> centers(mesh);

    _ComputeCenters<1, Dimension>::compute(centers, mesh);

    return centers;
}









/* TODO implement GS to compute volume and normal vector
template <unsigned int Dimension,typename IndexType, typename Real>
std::vector<Vertex<Dimension, Real>> GrammSchmidt(std::vector<Vertex<Dimension, Real>> vectors){

}
*/






template <unsigned int dim, unsigned int Dimension>
struct _ComputeMeasures{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>&,MeshElements<Dimension, IndexType, Real, Reserve...>&){
        static_assert (Dimension > 3,"The measure computation of mesh of dimension higher than 3 is not implemented yet.");
        throw std::runtime_error("The measure computation of mesh of dimension higher than 3 is not implemented yet.");
    }
};



template <>
struct _ComputeMeasures<3, 3>{
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

template <>
struct _ComputeMeasures<2, 2>{
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



template <>
struct _ComputeMeasures<2, 3>{
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







template <unsigned int Dimension>
struct _ComputeMeasures<1, Dimension>{
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





template <unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> ComputeMeasures(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
    MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>> measures(mesh);

    _ComputeMeasures<1, Dimension>::compute(measures, mesh);

    return measures;
}







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









template <unsigned int CurrentDimension, unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool End, bool Descend>
struct MeshRun {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){


        auto i = mesh.template getElements<CurrentDimension>().at(index);
        for (auto sube: mesh.template getElement<CurrentDimension>(i.getIndex()).getSubelements())
        MeshRun< CurrentDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == CurrentDimension - 1, Descend>::run(mesh, origElementIndex, sube.index, fun);


    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<MeshDimension, StartDimension, TargetDimension, MeshDimension, false, Descend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){

        auto& cell = mesh.getCells().at(index);
        IndexType tmpFace = cell.getBoundaryElementIndex();
        do {
            MeshRun<MeshDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == MeshDimension - 1, Descend>::run(mesh, origElementIndex, tmpFace, fun);
            tmpFace = mesh.getFaces().at(tmpFace).getNextBElem(cell.getIndex());
        } while (tmpFace != cell.getBoundaryElementIndex());

    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<1, StartDimension, TargetDimension, MeshDimension, false, Descend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){

        auto& edge = mesh.getEdges().at(index);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Descend>::run(mesh, origElementIndex, edge.getVertexAIndex(), fun);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Descend>::run(mesh, origElementIndex, edge.getVertexBIndex(), fun);
    }
};





template <unsigned int CurrentDimension,unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Descend>
struct MeshRun<CurrentDimension, StartDimension, TargetDimension, MeshDimension, true, Descend> {

    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& ,
                    IndexType origElementIndex,
                    IndexType index,
                    Functor fun){
        static_assert (std::is_assignable<std::function<void(IndexType, IndexType)>,Functor>::value,
                       "The Functor fun must be a function with void return type and two arguments of IndexType, the first is index of StartDimension element and the second is the index of the TargetDimension element");
        if(Descend){
            fun(origElementIndex, index);
        }else{
            fun(index, origElementIndex);
        }
    }
};





template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension>
struct MeshApply {
    template<typename Functor, typename IndexType, typename Real, unsigned int ...Reserve>
    static void apply(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                      Functor f) {
        for (auto& startElement : mesh.template getElements<(StartDimension > TargetDimension) ? StartDimension : TargetDimension>()){
            MeshRun<
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? TargetDimension : StartDimension,
                    MeshDimension,
                    StartDimension == TargetDimension,
                    (StartDimension > TargetDimension)>::run(mesh, startElement.getIndex(), startElement.getIndex(), f);
        }
    }
};






enum Order{
    ORDER_ASCEND,
    ORDER_ORIGINAL
};

template<unsigned int StartDim, unsigned int TargetDim, Order order = Order::ORDER_ASCEND>
struct MeshConnections {
    /**
     * @brief connections<HR>
     * Detects connections of mesh elements of StartDim to TargetDim.
     * Returns a MeshDataContainer of set<IndexType> allocated to StartDim elements.
     * The indexes are ordered in ascending way.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<std::set<IndexType>, StartDim> result(mesh);
        MeshApply<StartDim, TargetDim, MeshDimension>::apply(mesh, [&result](IndexType ori, IndexType element){
            result.template getDataByPos<0>().at(ori).insert(element);
        });

        return result;
    }
};


template<unsigned int StartDim, unsigned int TargetDim>
struct MeshConnections<StartDim, TargetDim, Order::ORDER_ORIGINAL> {

    /**
     * @brief orderedConnections<HR>
     * This function returns connection in original sequence as in the mesh.
     * @param mesh
     * @return
     */
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::vector<IndexType>, StartDim> connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<std::map<IndexType, IndexType>, StartDim> tempMap(mesh);

        MeshApply<StartDim, TargetDim, MeshDimension>::apply(mesh, [&tempMap](IndexType ori, IndexType element){
            IndexType size = tempMap.template getDataByPos<0>().at(ori).size();
            tempMap.template getDataByPos<0>().at(ori).insert({element, size});
        });

        MeshDataContainer<std::vector<IndexType>, StartDim> result(mesh);
        for (IndexType i = 0; i < mesh.template getElements<StartDim>().size(); i++){
            //resize the vector at the position
            result.template getDataByPos<0>().at(i).resize(
                tempMap.template getDataByPos<0>().at(i).size()
            );

            for(std::pair<const IndexType, IndexType>& mapElem : tempMap.template getDataByPos<0>().at(i)) {
                result.template getDataByPos<0>().at(i).at(mapElem.second) = mapElem.first;
            }
        }
        return result;
    }

};






template <unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim = StartDim, Order order = Order::ORDER_ASCEND>
class MeshNegborhood{
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> neighbors(
                MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::set<IndexType>, StartDim> result;
        MeshDataContainer<std::set<IndexType>, StartDim> firstConnections = MeshConnections<StartDim, ConnectingDim, order>::connections(mesh);
        MeshDataContainer<std::set<IndexType>, StartDim> secondConnections = MeshConnections<ConnectingDim, ConnectedDim, order>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<StartDim>().size(); elementIndex++) {

            for (IndexType& firstConectedElem : firstConnections.template getDataByPos<0>().at(elementIndex)){

                for (IndexType& neighborIndex : secondConnections.template getDataByPos<0>().at(firstConectedElem)){

                    if (StartDim == ConnectedDim && elementIndex == neighborIndex) {
                        continue;
                    } else {
                        result.template getDataByPos<0>().at(elementIndex).insert(neighborIndex);
                    }
                }
            }
        }

        return result;
    }

};


template <unsigned int StartDim, unsigned int ConnectingDim, unsigned int ConnectedDim>
class MeshNegborhood<StartDim, ConnectingDim, ConnectedDim, Order::ORDER_ORIGINAL>{
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> neighbors(
                MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {

        MeshDataContainer<std::map<IndexType, IndexType>, StartDim> tempResultMap(mesh);
        MeshDataContainer<std::vector<IndexType>, StartDim> firstConnections = MeshConnections<StartDim, ConnectingDim, Order::ORDER_ORIGINAL>::connections(mesh);
        MeshDataContainer<std::vector<IndexType>, StartDim> secondConnections = MeshConnections<ConnectingDim, ConnectedDim, Order::ORDER_ORIGINAL>::connections(mesh);

        for (IndexType elementIndex = 0; elementIndex < mesh.template getElements<StartDim>().size(); elementIndex++) {

            for (IndexType& firstConectedElem : firstConnections.template getDataByPos<0>().at(elementIndex)){

                for (IndexType& neighborIndex : secondConnections.template getDataByPos<0>().at(firstConectedElem)){

                    if (StartDim == ConnectedDim && elementIndex == neighborIndex) {
                        continue;
                    } else {
                        IndexType pos = tempResultMap.template getDataByPos<0>().at(elementIndex).size();
                        tempResultMap.template getDataByPos<0>().at(elementIndex).insert({neighborIndex, pos});
                    }
                }
            }
        }

        MeshDataContainer<std::vector<IndexType>, StartDim> result;

        for (IndexType i = 0; i < mesh.template getElements<StartDim>().size(); i++){
            //resize the vector at the position
            result.template getDataByPos<0>().at(i).resize(
                tempResultMap.template getDataByPos<0>().at(i).size()
            );

            for(std::pair<IndexType, IndexType>& mapElem : tempResultMap) {
                result.template getDataByPos<0>().at(i).at(mapElem.second) = mapElem.first;
            }
        }
        return result;
    }

};



template<unsigned int FromDim, unsigned int ToDim, bool Descend = true>
struct MeshColouring {

    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, FromDim> colour(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<unsigned int, FromDim> result(mesh);

        DBGMSG("starting the coloring procedure");
        unsigned int reserve = 16;
        MeshDataContainer<std::valarray<bool>, ToDim> attachedColours(mesh, std::valarray<bool>(false, reserve));



        for (auto& startElement : mesh.template getElements<FromDim>()){
            std::valarray<bool> possibleColours(true,reserve);
            MeshRun<FromDim, FromDim, ToDim, MeshDimension, false, true>::
                run(mesh,
                    startElement.getIndex(),
                    startElement.getIndex(),
                    [&possibleColours, &attachedColours](IndexType, IndexType element){
                        DBGTRY(possibleColours &= !attachedColours.template getDataByPos<0>().at(element);)
                    }
                );

            // Select the first possible colour
            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
                if (selectedColour == possibleColours.size()){
                    reserve *= 2;
                    DBGVAR(reserve);
                    for (std::valarray<bool>& attColour : attachedColours.template getDataByPos<0>()){
                        std::valarray<bool> newAttColour(false, reserve);
                        for (size_t i = 0; i < attColour.size(); i++){
                            newAttColour[i] = attColour[i];
                        }
                        attColour.swap(newAttColour);
                    }
                    break;
                }
            }
            result.template getDataByPos<0>().at(startElement.getIndex()) = selectedColour;
            MeshRun<FromDim, FromDim, ToDim, MeshDimension, false, true>::
                run(mesh,
                    startElement.getIndex(),
                    startElement.getIndex(),
                    [selectedColour, &attachedColours](IndexType, IndexType element){
                        DBGTRY(attachedColours.template getDataByPos<0>().at(element)[selectedColour] = true;)
                    }
                );
        }
        return result;
    }
};


template<unsigned int FromDim, unsigned int ToDim>
struct MeshColouring <FromDim, ToDim, false> {
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<unsigned int, FromDim> colour(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        // resulting container of colours
        MeshDataContainer<unsigned int, FromDim> result(mesh);

        DBGMSG("starting the coloring procedure");
        unsigned int reserve = 16;
        // allocates memory to the given dimension
        MeshDataContainer<std::valarray<bool>, ToDim> attachedColours(mesh, std::valarray<bool>(false, reserve));

        auto connections = MeshConnections<FromDim, ToDim>::connections(mesh);

        for (auto& startElement : mesh.template getElements<FromDim>()){
            std::valarray<bool> possibleColours(true,reserve);
            for (IndexType element : connections.at(startElement)){

                DBGTRY(possibleColours &= !attachedColours.template getDataByPos<0>().at(element);)

            }

            // Select the first possible colour
            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
                if (selectedColour == possibleColours.size()){
                    reserve *= 2;

                    // If the number of colours exceeds the number of
                    // allocated bits, then allocate twice as much memory
                    for (std::valarray<bool>& attColour : attachedColours.template getDataByPos<0>()){
                        std::valarray<bool> newAttColour(false, reserve);
                        for (size_t i = 0; i < attColour.size(); i++){
                            newAttColour[i] = attColour[i];
                        }
                        attColour.swap(newAttColour);
                    }
                    break;
                }
            }

            result.template getDataByPos<0>().at(startElement.getIndex()) = selectedColour;

            for (IndexType element : connections.at(startElement)){
                DBGTRY(attachedColours.template getDataByPos<0>().at(element)[selectedColour] = true;)
            }
        }
        return result;
    }
};


template <unsigned int FromDim, unsigned int ToDim>
struct ColourMesh{

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
static MeshDataContainer<unsigned int, FromDim> colour(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
        ){
    return MeshColouring<FromDim, ToDim, (FromDim > ToDim)>::colour(mesh);
}
};



template<typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<2, IndexType, Real, Reserve...>& mesh,
                typename MeshElements<2, IndexType, Real, Reserve...>::template ElementType<2>& face,
                typename MeshElements<2, IndexType, Real, Reserve...>::Edge& edge
                ) {

    Vertex<2, Real> AminC = mesh.getVertices().at(edge.getVertexAIndex()) - face.getCenter();
    Vertex<2, Real> BminC = mesh.getVertices().at(edge.getVertexBIndex()) - face.getCenter();



    Real res = AminC[0]*BminC[1]-BminC[0]*AminC[1];
    return res > 0;
    throw std::runtime_error("can not determine orientation of edge " +
                             std::to_string(edge.getIndex()) + " wrt face: " + std::to_string(face.getIndex()));
}

template<typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<3, IndexType, Real, Reserve...>& mesh,
                typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<2>& face,
                typename MeshElements<3, IndexType, Real, Reserve...>::Edge& edge,
                Vector<3, Real> faceNormal
                ) {

    Vertex<3, Real> AminC = mesh.getVertices().at(edge.getVertexAIndex()) - face.getCenter();
    Vertex<3, Real> BminC = mesh.getVertices().at(edge.getVertexBIndex()) - face.getCenter();

    Real res = Real(0);

    for (IndexType i = 0; i < 3; i++){
        IndexType ipo = (i+1)%(3);
        IndexType ipt = (i+2)%(3);
        res += AminC[i]*BminC[ipo]*faceNormal[ipt]-AminC[ipt]*BminC[ipo]*faceNormal[i];

    }
    return res > 0;
}

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh, IndexType faceIndex, IndexType edgeIndex) {

    typename MeshElements<MeshDimension, IndexType, Real, Reserve...>::Edge& edge = mesh.getEdges().at(edgeIndex);
    typename MeshElements<MeshDimension, IndexType, Real, Reserve...>::template ElementType<2>& face = mesh.template getElements<2>().at(faceIndex);

    return edgeIsLeft(mesh, face, edge);

}

template<typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<3, IndexType, Real, Reserve...>& mesh, IndexType faceIndex, IndexType edgeIndex) {

    typename MeshElements<3, IndexType, Real, Reserve...>::Edge& edge = mesh.getEdges().at(edgeIndex);
    typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<2>& face = mesh.template getElements<2>().at(faceIndex);

    auto normals = ComputeFaceNormals(mesh);

    return edgeIsLeft(mesh, face, edge, normals[face]);

}

template<typename IndexType, typename Real, unsigned int ...Reserve>
MeshDataContainer<std::vector<bool>, 2> edgesOrientation(MeshElements<3, IndexType, Real, Reserve...>& mesh) {

    MeshDataContainer<std::vector<bool>, 2> orientations(mesh);
    auto normals = ComputeFaceNormals(mesh);

    for (auto& face : mesh.getFaces()) {
        orientations[face].resize(face.getSubelements().getNumberOfSubElements());
        for (IndexType i = 0; i < face.getSubelements().getNumberOfSubElements(); i++){
            typename MeshElements<3, IndexType, Real, Reserve...>::Edge& edge = mesh.getEdges().at(face.getSubelements()[i].index);

            orientations[face][i] = edgeIsLeft(mesh, face, edge, normals[face]);
        }
    }
    return orientations;
}




#endif // MESH_FUNCTIONS_H
