#ifndef MESHDATACONTAINER_H
#define MESHDATACONTAINER_H

#include "MeshElement.h"
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
    static constexpr unsigned int dimensionIndex(){
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
        static void alocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());
            Alocator<pos - 1>::alocateMemory(parent, mesh);
        }

        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const DataType& initialValue) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);
            Alocator<pos - 1>::alocateMemory(parent, mesh, initialValue);
        }
    };

    template<typename dummy>
    struct Alocator<0, dummy>{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template getDataByPos<0>().resize(
                        mesh.template getElements<std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());

        }
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const DataType& initialValue) {

            parent.template getDataByPos<0>().resize(
                        mesh.template getElements<std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);

        }
    };




    /**
     * @brief data
     * A structure containing vectors of specified type
     * alocated to match the mesh elements.
     */
    _DataContainer<DataType, sizeof... (Dimensions) - 1> data;


public:

    template<unsigned int dim>
    std::vector<DataType>& getDataByDim(){
        return data._DataContainer<DataType, dimensionIndex<dim>()>::_data;
    }


    template<unsigned int pos>
    std::vector<DataType>& getDataByPos(){
        return data._DataContainer<DataType,pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& at(MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& operator[](MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }

    MeshDataContainer(){}

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        alocateData(mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh, std::integer_sequence<unsigned int,Dimensions...>, DataType){
        alocateData(mesh);
    }


    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataType& initialValue){
        alocateData(mesh, initialValue);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataType& initialValue,
                      std::integer_sequence<unsigned int,Dimensions...>,
                      DataType){
        alocateData(mesh, initialValue);
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void alocateData(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Alocator<sizeof... (Dimensions) - 1>::alocateMemory(*this, mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void alocateData(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh, const DataType& initialValue){
        Alocator<sizeof... (Dimensions) - 1>::alocateMemory(*this, mesh,initialValue);
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
    static constexpr unsigned int dimensionIndex(){
        return DimensionPos<dim, 0, std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>::res();
    }
public:

    template<unsigned int pos>
    using DataType = typename std::tuple_element<pos,std::tuple<DataTypes...>>::type;

    template<unsigned int Pos, typename Dummy = void>
    struct _DataContainer : _DataContainer<Pos - 1, Dummy>{
        std::vector<DataType<Pos>> _data;
    };

    template<typename Dummy>
    struct _DataContainer<0, Dummy>{
        std::vector<DataType<0>> _data;
    };

    template<unsigned int pos, typename _DataType, typename... _DataTypes>
    struct Alocator{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());
            Alocator<pos + 1, _DataTypes...>::alocateMemory(parent, mesh);
        }

        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const _DataType& initialValue,
                                  const _DataTypes&... values) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);
            Alocator<pos + 1, _DataTypes...>::alocateMemory(parent, mesh, values...);
        }
    };

    template<typename _DataType, typename... _DataTypes>
    struct Alocator<sizeof... (Dimensions) - 1, _DataType, _DataTypes...>{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template getDataByPos<sizeof... (Dimensions) - 1>().resize(
                        mesh.template getElements<std::get<sizeof... (Dimensions) - 1>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());

        }
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void alocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const _DataType& initialValue,
                                  const _DataTypes&...) {

            parent.template getDataByPos<sizeof... (Dimensions) - 1>().resize(
                        mesh.template getElements<std::get<sizeof... (Dimensions) - 1>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);

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
    std::vector<std::tuple_element_t<dimensionIndex<dim>(), std::tuple<DataTypes...>>>& getDataByDim(){
        return data._DataContainer<dimensionIndex<dim>()>::_data;
    }


    /**
     * @brief GetDataPos
     * @return
     */
    template<unsigned int pos>
    std::vector<DataType<pos>>& getDataByPos(){
        return data._DataContainer<pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>& at(MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>& operator[](MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }



    MeshDataContainer(){}

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        alocateData(mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataTypes&... initialValues){
        alocateData(mesh, initialValues...);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void alocateData(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Alocator<0, DataTypes...>::alocateMemory(*this, mesh);
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void alocateData(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                     const DataTypes&... initialValues){
        Alocator<0, DataTypes...>::alocateMemory(*this, mesh, initialValues...);
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