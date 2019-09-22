#ifndef MESHDATACONTAINER_H
#define MESHDATACONTAINER_H

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


#endif // MESHDATACONTAINER_H
