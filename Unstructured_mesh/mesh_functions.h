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



template <unsigned int dim, unsigned int Dimension, unsigned int... DataDimensions>
struct _ComputeCenters{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vertex<Dimension, Real>, DataDimensions...>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

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
        _ComputeCenters<dim + 1, Dimension, DataDimensions...>::compute(centers, mesh);
    }

};

template <unsigned int Dimension, unsigned int... DataDimensions>
struct _ComputeCenters<Dimension, Dimension, DataDimensions...>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vertex<Dimension, Real>, DataDimensions...>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

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

template <unsigned int Dimension, unsigned int... DataDimensions>
struct _ComputeCenters<1, Dimension, DataDimensions...>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    static void compute(MeshDataContainer<Vertex<Dimension, Real>, DataDimensions...>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        std::vector<Vertex<Dimension, Real>>& edgeCenters = centers.template GetDataDim<1>();

        for (auto& edge : mesh.template GetElements<1>()) {

            edgeCenters.at(edge.GetIndex()) = (mesh.template GetElements<0>().at(edge.GetVertexAIndex()) +
                                mesh.template GetElements<0>().at(edge.GetVertexBIndex())) * 0.5;
        }

        DBGMSG("1");
        _ComputeCenters<2, Dimension, DataDimensions...>::compute(centers, mesh);
    }
};




template <unsigned int Dimension,typename IndexType, typename Real, unsigned int ...Reserve>
auto __ComputeCenters(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
/*
    MeshDataContainer centers(mesh, std::make_integer_sequence<unsigned int, Dimension + 1>{},Vertex<Dimension, Real>{});

    _ComputeCenters<1, Dimension>::compute(centers, mesh);

    return centers;*/
    return 0;
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
