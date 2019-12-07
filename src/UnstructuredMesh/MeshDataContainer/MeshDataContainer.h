#ifndef MESHDATACONTAINER_H
#define MESHDATACONTAINER_H

#include "../MeshElements/MeshElement.h"



template<typename DataType, unsigned int Position, unsigned int MappedDimenion>
struct DataContainer : public std::vector<DataType> {
    using type = DataType;

    static constexpr unsigned int getPosition() {
        return Position;
    }

    static constexpr unsigned int getMappedDimension() {
        return MappedDimenion;
    }
};


/**
 * @brief The MeshDataContainer struct<HR>
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

    template<unsigned int pos>
    static constexpr unsigned int dimensionAt(){
        return std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...});
    }
public:

    template<typename _DataType, unsigned int _Dim>
    struct _DataContainer : _DataContainer<_DataType,_Dim - 1> {
        DataContainer<_DataType, _Dim, dimensionAt<_Dim>()> _data;
    };

    template<typename _DataType>
    struct _DataContainer<_DataType, 0>{
        DataContainer<_DataType, 0, dimensionAt<0U>()> _data;
    };

    /**
     * Data container type according to pos
     */
    template <unsigned int Pos>
    using DataContainerType = DataContainer<DataType, 0, dimensionAt<Pos>()>;

    /**
     * @brief size<HR>
     * Returns the number of vectors contained in the MeshDataContainer.
    * @return
    */
   static constexpr unsigned int size() {
       return sizeof... (Dimensions);
   }

private:
    template<unsigned int pos, typename dummy = void>
    struct Allocator{
        MeshDataContainer<DataType, Dimensions...>& parent;
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());
            Allocator<pos - 1>::allocateMemory(parent, mesh);
        }

        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const DataType& initialValue) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);
            Allocator<pos - 1>::allocateMemory(parent, mesh, initialValue);
        }


        /**
         * @brief allocateMemory
         * allocates memory according to another MashDataContainer of the same type
         * @param parent
         * @param meshDataContainer
         */
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<DataType, Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByPos<pos>().size());

            Allocator<pos - 1>::allocateMemory(parent, meshDataContainer);
        }

        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<DataType, Dimensions...>& meshDataContainer,
                                  const DataType& initialValue) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByPos<pos>().size(),
                        initialValue);
            Allocator<pos - 1>::allocateMemory(parent, meshDataContainer, initialValue);
        }
    };

    template<typename dummy>
    struct Allocator<0, dummy>{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template getDataByPos<0>().resize(
                        mesh.template getElements<std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());

        }
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const DataType& initialValue) {

            parent.template getDataByPos<0>().resize(
                        mesh.template getElements<std::get<0>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);

        }

        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<DataType, Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByPos<0>().size());
        }

        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<DataType, Dimensions...>& meshDataContainer,
                                  const DataType& initialValue) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByPos<0>().size(),
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
    DataContainer<DataType, dimensionIndex<dim>(), dim>& getDataByDim(){
        return data._DataContainer<DataType, dimensionIndex<dim>()>::_data;
    }

    template<unsigned int dim>
    const DataContainer<DataType, dimensionIndex<dim>(), dim>& getDataByDim() const {
        return data._DataContainer<DataType, dimensionIndex<dim>()>::_data;
    }


    template<unsigned int pos>
    DataContainer<DataType, pos, dimensionAt<pos>()>& getDataByPos(){
        return data._DataContainer<DataType,pos>::_data;
    }

    template<unsigned int pos>
    const DataContainer<DataType, pos, dimensionAt<pos>()>& getDataByPos() const {
        return data._DataContainer<DataType,pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    const DataType& at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) const {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    const DataType& operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) const {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }


    MeshDataContainer(){}

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        allocateData(mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      std::integer_sequence<unsigned int,Dimensions...>, DataType){
        allocateData(mesh);
    }


    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataType& initialValue){
        allocateData(mesh, initialValue);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataType& initialValue,
                      std::integer_sequence<unsigned int,Dimensions...>,
                      DataType){
        allocateData(mesh, initialValue);
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh, const DataType& initialValue){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, mesh,initialValue);
    }


    void allocateData(const MeshDataContainer<DataType, Dimensions...>& meshDataContainer){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, meshDataContainer);
    }

    void allocateData(const MeshDataContainer<DataType, Dimensions...>& meshDataContainer, const DataType& initialValue){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, meshDataContainer,initialValue);
    }

    MeshDataContainer<DataType, Dimensions...>& operator=(const MeshDataContainer<DataType, Dimensions...>& rhs) {
        this->data = rhs.data;
        return *this;
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

    template<unsigned int pos>
    static constexpr unsigned int dimensionAt(){
        return std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...});
    }
public:

    template<unsigned int pos>
    using DataType = typename std::tuple_element<pos,std::tuple<DataTypes...>>::type;

    template<unsigned int Pos, typename Dummy = void>
    struct _DataContainer : _DataContainer<Pos - 1, Dummy>{
        DataContainer<DataType<Pos>, Pos, dimensionAt<Pos>()> _data;
        //std::vector<DataType<Pos>> _data;
    };

    template<typename Dummy>
    struct _DataContainer<0, Dummy>{
        DataContainer<DataType<0>, 0, dimensionAt<0>()> _data;
        //std::vector<DataType<0>> _data;
    };

    /**
     * Data container type according to pos
     */
    template <unsigned int Pos>
    using DataContainerType = DataContainer<DataType<Pos>, 0, dimensionAt<Pos>()>;

    /**
     * @brief size<HR>
     * Returns the number of vectors contained in the MeshDataContainer.
    * @return
    */
   static constexpr unsigned int size() {
       return sizeof... (Dimensions);
   }

    template<unsigned int pos, typename _DataType, typename... _DataTypes>
    struct Allocator{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());
            Allocator<pos + 1, _DataTypes...>::allocateMemory(parent, mesh);
        }

        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const _DataType& initialValue,
                                  const _DataTypes&... values) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);
            Allocator<pos + 1, _DataTypes...>::allocateMemory(parent, mesh, values...);
        }

        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,
                                  const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByPos<pos>().size());
            Allocator<pos + 1, _DataTypes...>::allocateMemory(parent, meshDataContainer);
        }

        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer,
                                  const _DataType& initialValue,
                                  const _DataTypes&... values) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByPos<pos>().size(),
                        initialValue);
            Allocator<pos + 1, _DataTypes...>::allocateMemory(parent, meshDataContainer, values...);
        }
    };

    template<typename _DataType, typename... _DataTypes>
    struct Allocator<sizeof... (Dimensions) - 1, _DataType, _DataTypes...>{
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template getDataByPos<sizeof... (Dimensions) - 1>().resize(
                        mesh.template getElements<std::get<sizeof... (Dimensions) - 1>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size());

        }
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const _DataType& initialValue,
                                  const _DataTypes&...) {

            parent.template getDataByPos<sizeof... (Dimensions) - 1>().resize(
                        mesh.template getElements<std::get<sizeof... (Dimensions) - 1>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...})>().size(),
                        initialValue);

        }


        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,
                                  const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByPos<0>().size());

        }

        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer,
                                  const _DataType& initialValue,
                                  const _DataTypes&...) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByPos<0>().size(),
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
    DataContainer<std::tuple_element_t<dimensionIndex<dim>(), std::tuple<DataTypes...>>, dimensionIndex<dim>(), dim>&
        getDataByDim(){
        return data._DataContainer<dimensionIndex<dim>()>::_data;
    }


    /**
     * @brief GetDataPos
     * @return
     */
    template<unsigned int pos>
    DataContainer<DataType<pos>, pos, dimensionAt<pos>()>& getDataByPos(){
        return data._DataContainer<pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>& at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>& operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }



    MeshDataContainer(){}

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        allocateData(mesh);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataTypes&... initialValues){
        allocateData(mesh, initialValues...);
    }

    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Allocator<0, DataTypes...>::allocateMemory(*this, mesh);
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                     const DataTypes&... initialValues){
        Allocator<0, DataTypes...>::allocateMemory(*this, mesh, initialValues...);
    }


    void allocateData(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer){
        Allocator<0, DataTypes...>::allocateMemory(*this, meshDataContainer);
    }



    void allocateData(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer,
                     const DataTypes&... initialValues){
        Allocator<0, DataTypes...>::allocateMemory(*this, meshDataContainer, initialValues...);
    }

    MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& operator=(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& rhs) {
        this->data = rhs.data;
        return *this;
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
