#ifndef MESH_FUNCTIONS_H
#define MESH_FUNCTIONS_H
#include "mesh_element.h"


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

};


template <unsigned int dim, unsigned int Dimension>
struct _ComputeCenters{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    void compute(MeshDataContainer<Vertex<Dimension, Real>, std::make_index_sequence<Dimension + 1>{}>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& elemCenters = centers.template GetDataDim<dim>();
        auto& subElemCenters = centers.template GetDataDim<dim - 1>();


        for (IndexType i = 0; i < mesh.template GetElements<1>().size(); i++) {
            auto& element = mesh.template GetElements<dim>().at(i);

            Real subElemCnt = 0;
            for(IndexType sub : elemCenters.GetSubelements()){
                elemCenters.at(i) +=  subElemCenters.at(sub);
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
    void compute(MeshDataContainer<Vertex<Dimension, Real>, std::make_index_sequence<Dimension + 1>{}>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto& elemCenters = centers.template GetDataDim<Dimension>();
        auto& subElemCenters = centers.template GetDataDim<Dimension - 1>();


        for (IndexType i = 0; i < mesh.template GetElements<1>().size(); i++) {
            auto& element = mesh.template GetElements<Dimension>().at(i);

            Real subElemCnt = 0;

            IndexType tmpFaceIndex = element.GetBElemIndex();
            do {
                elemCenters.at(i) +=  subElemCenters.at(tmpFaceIndex);
                subElemCnt++;
            } while (tmpFaceIndex != element.GetBElemIndex());



            elemCenters.at(i) /= subElemCnt;
        }

    }

};

template <unsigned int Dimension>
struct _ComputeCenters<1, Dimension>{
    template <typename IndexType, typename Real, unsigned int ...Reserve>
    void compute(MeshDataContainer<Vertex<Dimension, Real>, std::make_index_sequence<Dimension + 1>{}>& centers,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

        auto & edgeCenters = centers.template GetDataDim<1>();

        for (IndexType i = 0; i < mesh.template GetElements<1>().size(); i++) {
            auto& edge = mesh.template GetElements<1>().at(i);
            edgeCenters.at(i) = 0.5 * (mesh.template GetElements<0>().at(edge.GetVertexAIndex) +
                                mesh.template GetElements<0>().at(edge.GetVertexAIndex));
        }

        _ComputeCenters<2, Dimension>::compute(centers, mesh);
    }

};

template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
MeshDataContainer<Vertex<Dimension, Real>, std::make_integer_sequence<unsigned int, Dimension + 1>{}> ComputeCenters(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

    MeshDataContainer<Vertex<Dimension, Real>, std::make_index_sequence<Dimension + 1>{}> centers(mesh);

    _ComputeCenters<1, Dimension>::compute(centers, mesh);

    return centers;
}

#endif // MESH_FUNCTIONS_H
