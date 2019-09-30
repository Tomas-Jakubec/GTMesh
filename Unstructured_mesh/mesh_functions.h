
#ifndef MESH_FUNCTIONS_H
#define MESH_FUNCTIONS_H
#include "mesh_element.h"
#include "meshdatacontainer.h"
#include "vector.h"
#include <valarray>
#include <set>

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

        auto& elemCenters = centers.template GetDataDim<dim>();
        auto& subElemCenters = centers.template GetDataDim<dim - 1>();


        for (IndexType i = 0; i < mesh.template getElements<dim>().size(); i++) {
            auto& element = mesh.template getElements<dim>().at(i);

            Real subElemCnt = 0;
            for(auto& sub : element.getSubelements()){
                elemCenters.at(i) +=  subElemCenters.at(sub.index);
                subElemCnt++;
            }

            elemCenters.at(i) /= subElemCnt;
        }

        DBGMSG(dim);
        _ComputeCenters<dim + 1, Dimension>::compute(centers, mesh);
    }
};

template <unsigned int Dimension>
struct _ComputeCenters<Dimension, Dimension>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& elemCenters = centers.template GetDataDim<Dimension>();
        auto& subElemCenters = centers.template GetDataDim<Dimension - 1>();


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
        DBGMSG(Dimension);
    }

};

template <unsigned int Dimension>
struct _ComputeCenters<1, Dimension>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(
            MakeMeshDataContainer_t<Vertex<Dimension, Real>, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& centers,
            MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        std::vector<Vertex<Dimension, Real>>& edgeCenters = centers.template GetDataDim<1>();

        for (auto& edge : mesh.template getElements<1>()) {

            edgeCenters.at(edge.getIndex()) = (mesh.template getElements<0>().at(edge.getVertexAIndex()) +
                                mesh.template getElements<0>().at(edge.getVertexBIndex())) * 0.5;
        }

        DBGMSG("1");
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

        auto& cellMeasures = measures.template GetDataDim<3>();

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

                Real tmp = distance * measures.template GetDataDim<2>().at(tmpFace);
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

        auto& surfaceMeasures = measures.template GetDataDim<2>();

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

        auto& surfaceMeasures = measures.template GetDataDim<2>();

        for (typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<2>& face : mesh.template getElements<2>()) {

            Real measure = Real();
            Vertex<3,Real>& faceCenter = face.getCenter();
            for(auto sube : face.getSubelements()){

                Vertex<3,Real>& a = mesh.getVertices().at(mesh.getEdges().at(sube.index).getVertexAIndex());
                Vertex<3,Real>& b = mesh.getVertices().at(mesh.getEdges().at(sube.index).getVertexBIndex());

                Real distance = Real();

                Real param = -1.0*(((a-faceCenter)*(b-a))/((b-a).sumOfSquares()));

                distance = (a-faceCenter+(b-a)*param).normEukleid();

                Real tmp = distance * measures.template GetDataDim<1>().at(sube.index);
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

        auto& edgeLengths = measures.template GetDataDim<1>();

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






namespace temp1 {

template <unsigned int CurrentDimension, unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool End, bool Ascend>
struct MeshRun {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){


        auto i = mesh.template getElements<CurrentDimension>().at(index);
        for (auto sube: mesh.template getElement<CurrentDimension>(i.getIndex()).getSubelements())
        MeshRun< CurrentDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == CurrentDimension - 1, Ascend>::run(mesh, origElementIndex, sube.index, fun);


    }
};

template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Ascend>
struct MeshRun<MeshDimension, StartDimension, TargetDimension, MeshDimension, false, Ascend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){

        auto& cell = mesh.getCells().at(index);
        IndexType tmpFace = cell.getBoundaryElementIndex();
        do {
            MeshRun<MeshDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == MeshDimension - 1, Ascend>::run(mesh, origElementIndex, tmpFace, fun);
            tmpFace = mesh.getFaces().at(tmpFace).getNextBElem(cell.getIndex());
        } while (tmpFace != cell.getBoundaryElementIndex());

    }
};


template <unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Ascend>
struct MeshRun<1, StartDimension, TargetDimension, MeshDimension, false, Ascend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){

        auto& edge = mesh.getEdges().at(index);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Ascend>::run(mesh, origElementIndex, edge.getVertexAIndex(), fun);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Ascend>::run(mesh, origElementIndex, edge.getVertexBIndex(), fun);
    }
};



template <unsigned int CurrentDimension,unsigned int StartDimension, unsigned int TargetDimension, unsigned int MeshDimension, bool Ascend>
struct MeshRun<CurrentDimension, StartDimension, TargetDimension, MeshDimension, true, Ascend> {

    template<typename Func, typename IndexType, typename Real, unsigned int ...Reserve>
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& ,
                    IndexType origElementIndex,
                    IndexType index,
                    Func fun){
        if(Ascend){
            fun(StartDimension, TargetDimension, origElementIndex, index);
        }else{
            fun(TargetDimension, StartDimension, index, origElementIndex);
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
                    false,
                    (StartDimension > TargetDimension)>::run(mesh, startElement.getIndex(), startElement.getIndex(), f);
        }
    }
};


template<unsigned int StartDim, unsigned int TargetDim>
struct MeshConnections {
    template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshDataContainer<std::set<IndexType>, StartDim> connections(
            MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh
            ) {
        MeshDataContainer<std::set<IndexType>, StartDim> result(mesh);
        MeshApply<StartDim, TargetDim, MeshDimension>::apply(mesh, [&result](unsigned int, unsigned int, IndexType ori, IndexType element){
            result.template GetDataPos<0>().at(ori).insert(element);
        });

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
                    [&possibleColours, &attachedColours](unsigned int, unsigned int, IndexType, IndexType element){
                        DBGTRY(possibleColours &= !attachedColours.template GetDataPos<0>().at(element);)
                    }
                );

            // Select the first possible colour
            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
            }
            result.template GetDataPos<0>().at(startElement.getIndex()) = selectedColour;
            MeshRun<FromDim, FromDim, ToDim, MeshDimension, false, true>::
                run(mesh,
                    startElement.getIndex(),
                    startElement.getIndex(),
                    [selectedColour, &attachedColours](unsigned int, unsigned int, IndexType, IndexType element){
                        DBGTRY(attachedColours.template GetDataPos<0>().at(element)[selectedColour] = true;)
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
        MeshDataContainer<unsigned int, FromDim> result(mesh);

        DBGMSG("starting the coloring procedure");
        unsigned int reserve = 16;
        MeshDataContainer<std::valarray<bool>, ToDim> attachedColours(mesh, std::valarray<bool>(false, reserve));

        auto connections = MeshConnections<FromDim, ToDim>::connections(mesh);

        for (auto& startElement : mesh.template getElements<FromDim>()){
            std::valarray<bool> possibleColours(true,reserve);
            for (IndexType element : connections.at(startElement)){

                DBGTRY(possibleColours &= !attachedColours.template GetDataPos<0>().at(element);)

            }

            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
            }

            result.template GetDataPos<0>().at(startElement.getIndex()) = selectedColour;

            for (IndexType element : connections.at(startElement)){
                DBGTRY(attachedColours.template GetDataPos<0>().at(element)[selectedColour] = true;)
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

}


#endif // MESH_FUNCTIONS_H
