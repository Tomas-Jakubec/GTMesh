#ifndef MESH_FUNCTIONS_H
#define MESH_FUNCTIONS_H
#include "mesh_element.h"
#include "../debug/debug.h"


/**
 * @brief The MeshDataContainer struct
 *
 * A struct designed to manage data boud to mesh.
 * Creates a serie of vectors sized acording to dimension.
 */
template <typename DataType, unsigned int ...Dimensions>
struct MeshDataContainer{
private:

    template<unsigned int dim, unsigned int pos, unsigned int _dim>
    struct DimensionPos : DimensionPos<dim, pos + 1,std::get<pos + 1>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>{};

    template<unsigned int dim, unsigned int pos>
    struct DimensionPos<dim, pos, dim>{
        static constexpr unsigned int res(){return pos;}
    };


    template<unsigned int dim>
    static constexpr unsigned int DimIndex(){
        return DimensionPos<dim, 0, std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>::res();
    }

    template<typename _DataType, unsigned int _Dim>
    struct _DataContainer : _DataContainer<_DataType,_Dim - 1>{
        std::vector<_DataType> _data;
    };

    template<typename _DataType>
    struct _DataContainer<_DataType, 0>{
        std::vector<_DataType> _data;
    };

    template<unsigned int pos, typename dummy = void>
    struct Alocator{
        MeshDataContainer<DataType, Dimensions...>& parent;
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void AlocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template GetDataPos<pos>().resize(
                        mesh.template GetElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());
            Alocator<pos - 1>::AlocateMemory(parent, mesh);
        }
    };

    template<typename dummy>
    struct Alocator<0, dummy>{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void AlocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template GetDataPos<0>().resize(
                        mesh.template GetElements<std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());

        }
    };




    /**
     * @brief data
     * A structure containing vectors of specified type
     * alocated to match the mesh elements.
     */
    _DataContainer<DataType, sizeof... (Dimensions) - 1> data;


public:

    /**
     * @brief GetDataDim
     * @return
     */
    template<unsigned int dim>
    std::vector<DataType>& GetDataDim(){
        return data._DataContainer<DataType, DimIndex<dim>()>::_data;
    }


    /**
     * @brief GetDataPos
     * @return
     */
    template<unsigned int pos>
    std::vector<DataType>& GetDataPos(){
        return data._DataContainer<DataType,pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& at(MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return GetDataDim<ElementDim>().at(element.GetIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& operator[](MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return GetDataDim<ElementDim>()[element.GetIndex()];
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Alocator<sizeof... (Dimensions) - 1>::AlocateMemory(*this, mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh, std::integer_sequence<unsigned int,Dimensions...>, DataType){
        DBGVAR(sizeof... (Dimensions))
        Alocator<sizeof... (Dimensions) - 1>::AlocateMemory(*this, mesh);
    }
};




/**
 * @brief The MeshDataContainer struct
 *
 * A struct designed to manage data boud to mesh.
 * Creates a serie of vectors sized acording to dimension.
 */
template <typename ...DataTypes, unsigned int ...Dimensions>
struct MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>{
private:

    template<unsigned int dim, unsigned int pos, unsigned int _dim>
    struct DimensionPos : DimensionPos<dim, pos + 1,std::get<pos + 1>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>{};

    template<unsigned int dim, unsigned int pos>
    struct DimensionPos<dim, pos, dim>{
        static constexpr unsigned int res(){return pos;}
    };

public:
    template<unsigned int dim>
    static constexpr unsigned int DimIndex(){
        return DimensionPos<dim, 0, std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>::res();
    }
public:

    template<unsigned int pos>
    using DataType = typename std::tuple_element<pos,std::tuple<DataTypes...>>::type;

    template<unsigned int _Dim, typename Dummy = void>
    struct _DataContainer : _DataContainer<_Dim - 1, Dummy>{
        std::vector<DataType<_Dim>> _data;
    };

    template<typename Dummy>
    struct _DataContainer<0, Dummy>{
        std::vector<DataType<0>> _data;
    };

    template<unsigned int pos, typename dummy = void>
    struct Alocator{
        MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent;
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void AlocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template GetDataPos<pos>().resize(
                        mesh.template GetElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());
            Alocator<pos - 1>::AlocateMemory(parent, mesh);
        }
    };

    template<typename dummy>
    struct Alocator<0, dummy>{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void AlocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template GetDataPos<0>().resize(
                        mesh.template GetElements<std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());

        }
    };




    /**
     * @brief data
     * A structure containing vectors of specified type
     * alocated to match the mesh elements.
     */
    _DataContainer<sizeof... (Dimensions) - 1> data;


public:

    /**
     * @brief GetDataDim
     * @return
     */
    template<unsigned int dim>
    std::vector<std::tuple_element_t<DimIndex<dim>(), std::tuple<DataTypes...>>>& GetDataDim(){
        return data._DataContainer<DimIndex<dim>()>::_data;
    }


    /**
     * @brief GetDataPos
     * @return
     */
    template<unsigned int pos>
    std::vector<DataType<pos>>& GetDataPos(){
        return data._DataContainer<pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<DimIndex<ElementDim>(), std::tuple<DataTypes...>>& at(MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return GetDataDim<ElementDim>().at(element.GetIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<DimIndex<ElementDim>(), std::tuple<DataTypes...>>& operator[](MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return GetDataDim<ElementDim>()[element.GetIndex()];
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Alocator<sizeof... (Dimensions) - 1>::AlocateMemory(*this, mesh);
    }

};




/**
 * MakeMeshDataContainer
 */
template<typename... Params>
struct MakeMeshDataContainer {};

template<typename Type, unsigned int... Dimensions>
struct MakeMeshDataContainer<Type, std::integer_sequence<unsigned int, Dimensions...>>{
    using type = MeshDataContainer<Type, Dimensions...>;
};


template<typename... Types, unsigned int... Dimensions>
struct MakeMeshDataContainer<std::tuple<Types...>, std::integer_sequence<unsigned int, Dimensions...>>{
    using type = MeshDataContainer<std::tuple<Types...>, Dimensions...>;
};


template<typename... T>
using MakeMeshDataContainer_t = typename MakeMeshDataContainer<T...>::type;


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
            for(auto& sub : element.GetSubelements()){
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
            IndexType tmpFaceIndex = element.GetBoundaryElementIndex();
            do {
                elemCenters.at(i) +=  subElemCenters.at(tmpFaceIndex);
                subElemCnt++;
                tmpFaceIndex = mesh.GetFaces()[tmpFaceIndex].GetNextBElem(i);
            } while (tmpFaceIndex != element.GetBoundaryElementIndex());

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
            IndexType tmpFace = cell.GetBoundaryElementIndex();
            Real measure = Real();
            Vertex<3,Real>& cellCenter = cell.GetCenter();

            do {
                // select 3 different vertices
                IndexType vAIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).GetSubelements()[0].index).GetVertexAIndex();
                IndexType vBIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).GetSubelements()[0].index).GetVertexBIndex();
                IndexType vCIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).GetSubelements()[1].index).GetVertexAIndex();
                if(vCIndex == vAIndex || vCIndex == vBIndex) {
                    vCIndex = mesh.GetEdges().at(mesh.GetFaces().at(tmpFace).GetSubelements()[1].index).GetVertexBIndex();
                }

                Vertex<3,Real>& a = mesh.GetVertices().at(vAIndex);
                Vertex<3,Real>& b = mesh.GetVertices().at(vBIndex);
                Vertex<3,Real>& c = mesh.GetVertices().at(vCIndex);

                // preparing quiantities
                Vertex<3,Real> vAmcC = (a-cellCenter);
                Vertex<3,Real> vBmA = (b-a);
                Vertex<3,Real> vCmA = (c-a);
                Real inv_sqrBmA = 1.0 / vBmA.SumOfSquares();
                Real inv_sqrCmA = 1.0 / vCmA.SumOfSquares();

                Real denominator = 1.0 / (1.0 - (pow(vCmA*vBmA,2) * inv_sqrBmA * inv_sqrCmA));


                Real param_t = -denominator * (((vAmcC*vBmA) * inv_sqrBmA) - (inv_sqrBmA*inv_sqrCmA*(vAmcC * vCmA)*(vCmA*vBmA)));
                //param_t *= inv_sqrBmA;
                Real param_s = -denominator * (((vAmcC*vCmA) * inv_sqrCmA) - (inv_sqrBmA*inv_sqrCmA*(vAmcC * vBmA)*(vCmA*vBmA)));

                Real distance = (vAmcC + (vBmA * param_t) + (vCmA * param_s)).NormEukleid();

                Real tmp = distance * measures.template GetDataDim<2>().at(tmpFace);
                measure += tmp / 3.0;

                tmpFace = mesh.GetFaces().at(tmpFace).GetNextBElem(cell.GetIndex());
            } while (tmpFace != cell.GetBoundaryElementIndex());

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
            IndexType tmpEdge = cell.GetBoundaryElementIndex();
            Real measure = Real();
            Vertex<2,Real>& cellCenter = cell.GetCenter();
            do {
                Vertex<2,Real>& a = mesh.GetVertices().at(mesh.GetEdges().at(tmpEdge).GetVertexAIndex());
                Vertex<2,Real>& b = mesh.GetVertices().at(mesh.GetEdges().at(tmpEdge).GetVertexBIndex());
                double tmp = (cellCenter[0] - a[0]) * (b[1] - a[1]);
                tmp -= (cellCenter[1] - a[1]) * (b[0] - a[0]);
                measure += 0.5 * fabs(tmp);

                tmpEdge = mesh.GetEdges().at(tmpEdge).GetNextBElem(cell.GetIndex());
            } while (tmpEdge != cell.GetBoundaryElementIndex());

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
            Vertex<3,Real>& faceCenter = face.GetCenter();
            for(auto sube : face.GetSubelements()){

                Vertex<3,Real>& a = mesh.GetVertices().at(mesh.GetEdges().at(sube.index).GetVertexAIndex());
                Vertex<3,Real>& b = mesh.GetVertices().at(mesh.GetEdges().at(sube.index).GetVertexBIndex());

                Real distance = Real();

                Real param = -1.0*(((a-faceCenter)*(b-a))/((b-a).SumOfSquares()));

                distance = (a-faceCenter+(b-a)*param).NormEukleid();

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
                                               mesh.GetVertices().at(edge.GetVertexBIndex())).NormEukleid();
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















template<unsigned int MeshDimension, unsigned int ElementDim,typename IndexType, typename Real, unsigned int ...Reserve>
struct CellsVertices {
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh, IndexType index){
        DBGMSG("face number "<<index);
        for(auto i : mesh.template GetElement<ElementDim>(index).GetSubelements()) {
            CellsVertices<MeshDimension, ElementDim - 1, IndexType, Real, Reserve...>::run(mesh, i.index);
        }
    }
};


template<unsigned int MeshDimension,typename IndexType, typename Real, unsigned int ...Reserve>
struct CellsVertices<MeshDimension, MeshDimension, IndexType, Real, Reserve...> {
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){
        for(IndexType i = 0; i < mesh.GetCells().size(); i++){
            DBGMSG("cell number "<<i);
            for(auto j : mesh.template GetElement<MeshDimension>(i).GetSubelements()) {
                CellsVertices<MeshDimension, MeshDimension - 1, IndexType, Real, Reserve...>::run(mesh, j);
            }
        }
    }
};


template<unsigned int MeshDimension,typename IndexType, typename Real, unsigned int ...Reserve>
struct CellsVertices<MeshDimension, 1, IndexType, Real, Reserve...> {
    static void run(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh, IndexType index){

            auto e = mesh.template GetElement<1>(index);
            DBGVAR(mesh.GetVertices()[e.GetElement().GetVertexAIndex()], mesh.GetVertices()[e.GetElement().GetVertexBIndex()])

    }
};

#endif // MESH_FUNCTIONS_H
