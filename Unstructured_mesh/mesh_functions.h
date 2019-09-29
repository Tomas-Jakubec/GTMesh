
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


        for (IndexType i = 0; i < mesh.template GetElements<dim>().size(); i++) {
            auto& element = mesh.template GetElements<dim>().at(i);

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


        for (IndexType i = 0; i < mesh.template GetElements<Dimension>().size(); i++) {
            auto& element = mesh.template GetElements<Dimension>().at(i);

            Real subElemCnt = 0;
            IndexType tmpFaceIndex = element.getBoundaryElementIndex();
            do {
                elemCenters.at(i) +=  subElemCenters.at(tmpFaceIndex);
                subElemCnt++;
                tmpFaceIndex = mesh.GetFaces()[tmpFaceIndex].getNextBElem(i);
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

        for (auto& edge : mesh.template GetElements<1>()) {

            edgeCenters.at(edge.GetIndex()) = (mesh.template GetElements<0>().at(edge.GetVertexAIndex()) +
                                mesh.template GetElements<0>().at(edge.GetVertexBIndex())) * 0.5;
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

        for (typename MeshElements<3, IndexType, Real, Reserve...>::template ElemType<3>::type& cell : mesh.GetCells()) {
            IndexType tmpFace = cell.getBoundaryElementIndex();
            Real measure = Real();
            Vertex<3,Real>& cellCenter = cell.getCenter();

            do {
                // select 3 different vertices
                IndexType vAIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).getSubelements()[0].index).GetVertexAIndex();
                IndexType vBIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).getSubelements()[0].index).GetVertexBIndex();
                IndexType vCIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).getSubelements()[1].index).GetVertexAIndex();
                if(vCIndex == vAIndex || vCIndex == vBIndex) {
                    vCIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).getSubelements()[1].index).GetVertexBIndex();
                }

                Vertex<3,Real>& a = mesh.GetVertices().at(vAIndex);
                Vertex<3,Real>& b = mesh.GetVertices().at(vBIndex);
                Vertex<3,Real>& c = mesh.GetVertices().at(vCIndex);

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

                tmpFace = mesh.GetFaces().at(tmpFace).getNextBElem(cell.GetIndex());
            } while (tmpFace != cell.getBoundaryElementIndex());

            cellMeasures.at(cell.GetIndex()) = measure;
        }
    }
};

template <>
struct _ComputeMeasures<2, 2>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 2>>& measures,MeshElements<2, IndexType, Real, Reserve...>& mesh){

        auto& surfaceMeasures = measures.template GetDataDim<2>();

        for (typename MeshElements<2, IndexType, Real, Reserve...>::template ElemType<2>::type& cell : mesh.GetCells()) {
            IndexType tmpEdge = cell.getBoundaryElementIndex();
            Real measure = Real();
            Vertex<2,Real>& cellCenter = cell.getCenter();
            do {
                Vertex<2,Real>& a = mesh.GetVertices().at(mesh.GetEdges().at(tmpEdge).GetVertexAIndex());
                Vertex<2,Real>& b = mesh.GetVertices().at(mesh.GetEdges().at(tmpEdge).GetVertexBIndex());
                double tmp = (cellCenter[0] - a[0]) * (b[1] - a[1]);
                tmp -= (cellCenter[1] - a[1]) * (b[0] - a[0]);
                measure += 0.5 * fabs(tmp);

                tmpEdge = mesh.GetEdges().at(tmpEdge).getNextBElem(cell.GetIndex());
            } while (tmpEdge != cell.getBoundaryElementIndex());

            surfaceMeasures.at(cell.GetIndex()) = measure;
        }
    }
};



template <>
struct _ComputeMeasures<2, 3>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, 3>>& measures,MeshElements<3, IndexType, Real, Reserve...>& mesh){

        auto& surfaceMeasures = measures.template GetDataDim<2>();

        for (typename MeshElements<3, IndexType, Real, Reserve...>::template ElemType<2>::type& face : mesh.template GetElements<2>()) {

            Real measure = Real();
            Vertex<3,Real>& faceCenter = face.getCenter();
            for(auto sube : face.getSubelements()){

                Vertex<3,Real>& a = mesh.GetVertices().at(mesh.GetEdges().at(sube.index).GetVertexAIndex());
                Vertex<3,Real>& b = mesh.GetVertices().at(mesh.GetEdges().at(sube.index).GetVertexBIndex());

                Real distance = Real();

                Real param = -1.0*(((a-faceCenter)*(b-a))/((b-a).sumOfSquares()));

                distance = (a-faceCenter+(b-a)*param).normEukleid();

                Real tmp = distance * measures.template GetDataDim<1>().at(sube.index);
                measure += tmp * 0.5;
            }
            surfaceMeasures.at(face.GetIndex()) = measure;
        }
        _ComputeMeasures<3, 3>::compute(measures, mesh);
    }
};







template <unsigned int Dimension>
struct _ComputeMeasures<1, Dimension>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MakeMeshDataContainer_t<Real, make_custom_integer_sequence_t<unsigned int, 1, Dimension>>& measures,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& edgeLengths = measures.template GetDataDim<1>();

        for (auto& edge : mesh.GetEdges()) {
            edgeLengths.at(edge.GetIndex()) = (mesh.GetVertices().at(edge.GetVertexAIndex()) -
                                               mesh.GetVertices().at(edge.GetVertexBIndex())).normEukleid();
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
        for (auto& face : mesh.GetEdges()) {
            Vertex<2,Real> a = mesh.GetVertices().at(face.GetVertexAIndex());
            Vertex<2,Real> b = mesh.GetVertices().at(face.GetVertexBIndex());
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
        for (auto& face : mesh.GetFaces()) {

            bool vectorSign = true;
            IndexType cellIndex = face.getCellLeftIndex();
            if (
                    cellIndex == INVALID_INDEX(IndexType) ||
                    (cellIndex & BOUNDARY_INDEX(IndexType)) == BOUNDARY_INDEX(IndexType)
                ) {
                vectorSign = false;
                cellIndex = face.getCellRightIndex();
            }

            Vertex<3,Real>& cellCenter = mesh.GetCells().at(cellIndex).getCenter();


            // select 3 different vertices
            IndexType vAIndex = mesh.GetEdges().at(face.getSubelements()[0].index).GetVertexAIndex();
            IndexType vBIndex = mesh.GetEdges().at(face.getSubelements()[0].index).GetVertexBIndex();
            IndexType vCIndex = mesh.GetEdges().at(face.getSubelements()[1].index).GetVertexAIndex();
            if(vCIndex == vAIndex || vCIndex == vBIndex) {
                vCIndex = mesh.GetEdges().at(face.getSubelements()[1].index).GetVertexBIndex();
            }

            Vertex<3,Real>& a = mesh.GetVertices().at(vAIndex);
            Vertex<3,Real>& b = mesh.GetVertices().at(vBIndex);
            Vertex<3,Real>& c = mesh.GetVertices().at(vCIndex);

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

    if(mesh.GetBoundaryCells().empty()){
        for(auto& face : mesh.GetFaces()) {
            if (face.getCellLeftIndex() != INVALID_INDEX(IndexType) &&
                face.getCellRightIndex() != INVALID_INDEX(IndexType)){

                distances.at(face) = (mesh.GetCells().at(face.getCellLeftIndex()).getCenter() -
                                      mesh.GetCells().at(face.getCellRightIndex()).getCenter()).normEukleid();

            } else if(face.getCellLeftIndex() != INVALID_INDEX(IndexType) &&
                      face.getCellRightIndex() == INVALID_INDEX(IndexType)){

                distances.at(face) = (mesh.GetCells().at(face.getCellLeftIndex()).getCenter() -
                                      face.getCenter()).normEukleid();

            } else if(face.getCellLeftIndex() == INVALID_INDEX(IndexType) &&
                      face.getCellRightIndex() != INVALID_INDEX(IndexType)){

                distances.at(face) = (mesh.GetCells().at(face.getCellRightIndex()).getCenter() -
                                      face.getCenter()).normEukleid();
            }
        }

    } else {

        for(auto& face : mesh.GetFaces()) {
            auto& cellLeft = (face.getCellLeftIndex() & BOUNDARY_INDEX(IndexType)) == BOUNDARY_INDEX(IndexType)?
                      mesh.GetBoundaryCells().at(face.getCellLeftIndex()&EXTRACTING_INDEX(IndexType)):
                      mesh.GetCells().at(face.getCellLeftIndex());
            auto& cellRight = (face.getCellRightIndex() & BOUNDARY_INDEX(IndexType)) == BOUNDARY_INDEX(IndexType)?
                      mesh.GetBoundaryCells().at(face.getCellRightIndex()&EXTRACTING_INDEX(IndexType)):
                      mesh.GetCells().at(face.getCellRightIndex());

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


        auto i = mesh.template GetElements<CurrentDimension>().at(index);
        for (auto sube: mesh.template GetElement<CurrentDimension>(i.GetIndex()).getSubelements())
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

        auto& cell = mesh.GetCells().at(index);
        IndexType tmpFace = cell.getBoundaryElementIndex();
        do {
            MeshRun<MeshDimension - 1, StartDimension, TargetDimension, MeshDimension, TargetDimension == MeshDimension - 1, Ascend>::run(mesh, origElementIndex, tmpFace, fun);
            tmpFace = mesh.GetFaces().at(tmpFace).getNextBElem(cell.GetIndex());
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

        auto& edge = mesh.GetEdges().at(index);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Ascend>::run(mesh, origElementIndex, edge.GetVertexAIndex(), fun);
        MeshRun<0, StartDimension, TargetDimension, MeshDimension, TargetDimension == 0, Ascend>::run(mesh, origElementIndex, edge.GetVertexBIndex(), fun);
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
        for (auto& startElement : mesh.template GetElements<(StartDimension > TargetDimension) ? StartDimension : TargetDimension>()){
            MeshRun<
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? StartDimension : TargetDimension,
                    (StartDimension > TargetDimension) ? TargetDimension : StartDimension,
                    MeshDimension,
                    false,
                    (StartDimension > TargetDimension)>::run(mesh, startElement.GetIndex(), startElement.GetIndex(), f);
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



        for (auto& startElement : mesh.template GetElements<FromDim>()){
            std::valarray<bool> possibleColours(true,reserve);
            MeshRun<FromDim, FromDim, ToDim, MeshDimension, false, true>::
                run(mesh,
                    startElement.GetIndex(),
                    startElement.GetIndex(),
                    [&possibleColours, &attachedColours](unsigned int, unsigned int, IndexType, IndexType element){
                        DBGTRY(possibleColours &= !attachedColours.template GetDataPos<0>().at(element);)
                    }
                );

            // Select the first possible colour
            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
            }
            result.template GetDataPos<0>().at(startElement.GetIndex()) = selectedColour;
            MeshRun<FromDim, FromDim, ToDim, MeshDimension, false, true>::
                run(mesh,
                    startElement.GetIndex(),
                    startElement.GetIndex(),
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

        for (auto& startElement : mesh.template GetElements<FromDim>()){
            std::valarray<bool> possibleColours(true,reserve);
            for (IndexType element : connections.at(startElement)){

                DBGTRY(possibleColours &= !attachedColours.template GetDataPos<0>().at(element);)

            }

            unsigned int selectedColour = 0;
            while (!possibleColours[selectedColour]) {
                selectedColour++;
            }

            result.template GetDataPos<0>().at(startElement.GetIndex()) = selectedColour;

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
