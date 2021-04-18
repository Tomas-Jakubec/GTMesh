#ifndef MESHDATACONTAINER_H
#define MESHDATACONTAINER_H

#include "../MeshElements/MeshElements.h"

// TODO move
/**
 * @brief The DimensionPos struct
 * realizes the method indexof in the parameter pack @a Dimensions.
 * If the searched value is nor present, the index @a pos ran out of bounds.
 */
template<unsigned int dim, unsigned int pos, unsigned int _dim, unsigned int... Dimensions>
struct DimensionPos : public std::conditional_t<dim == _dim,
                                                std::integral_constant<unsigned int, pos>,
                                                DimensionPos<dim, pos + 1, Dimensions...>>
{};

template<unsigned int dim, unsigned int pos, unsigned int _dim>
struct DimensionPos<dim, pos, _dim> : public std::integral_constant<unsigned int, pos>{
    static_assert(dim == _dim, "The desired number was not found.");
};

template<typename DataType, unsigned int MappedDimenion = 0>
struct DataContainer : public std::vector<DataType> {
    using type = DataType;


    static constexpr unsigned int getMappedDimension() {
        return MappedDimenion;
    }
};


/**
 * @brief The MeshDataContainer struct<HR>
 * A struct designed to manage data bound to a mesh.
 * Creates a serie of vectors sized acording to dimension.
 */
template <typename DataType, unsigned int ...Dimensions>
struct MeshDataContainer{
public:
    template<unsigned int dim>
    static constexpr unsigned int dimensionIndex(){
        return DimensionPos<dim, 0, Dimensions...>::value;
    }
    template<unsigned int pos>
    static constexpr unsigned int dimensionAt(){
        return std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...});
    }
private:

    /**
     * @brief The _DataContainer struct
     * constructs a system of DataContainers of given Dimension indexes.
     */
    template<typename _DataType, unsigned int Pos>
    struct _DataContainer : _DataContainer<_DataType,Pos - 1> {
        DataContainer<_DataType, dimensionAt<Pos>()> _data;
    };

    /**
     * @brief The _DataContainer<_DataType, _Tp2> struct
     * cycle terminating specialization
     */
    template<typename _DataType>
    struct _DataContainer<_DataType, 0>{
        DataContainer<_DataType, dimensionAt<0U>()> _data;
    };
public:
    /**
     * Data container type according to pos
     */
    template <unsigned int Pos>
    using DataContainerType = DataContainer<DataType, dimensionAt<Pos>()>;

    /**
     * @brief Returns the number of vectors contained in the MeshDataContainer.
    * @return
    */
   static constexpr unsigned int size() {
       return sizeof... (Dimensions);
   }

private:
   /**
    * @brief This struct recursively allocates data of MeshDataContainer.
    */
    template<unsigned int pos, typename dummy = void>
    struct Allocator{

        /**
         * @brief Resizes vector on positions according to the respective dimension
         * of mesh and initializes with defaul value if possible.
         * @param parent
         * @param mesh
         */
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<parent.template dimensionAt<pos>()>().size());
            Allocator<pos - 1>::allocateMemory(parent, mesh);
        }



        /**
         * @brief Resizes vector on positions according to the respective dimension
         * of mesh and initializes with @a initialValue.
         * @param parent
         * @param mesh
         * @param initialValue a value to be data initialized with
         */
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const DataType& initialValue) {
            parent.template getDataByPos<pos>().resize(
                        mesh.template getElements<parent.template dimensionAt<pos>()>().size(),
                        initialValue);
            Allocator<pos - 1>::allocateMemory(parent, mesh, initialValue);
        }


        /**
         * @brief Allocates memory according to another MashDataContainer which must
         * have data allocated to the @a Dimensions. The data is initialized by
         * default value if possible.
         * @param parent
         * @param meshDataContainer
         */
        template<typename _DataType, unsigned int ... _Dimensions>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<_DataType, _Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<pos>()>().size());

            Allocator<pos - 1>::allocateMemory(parent, meshDataContainer);
        }



        /**
         * @brief Allocates memory according to another %MashDataContainer which must
         * have data allocated to the @a Dimensions. The data is initialized by
         * @a initialValue.
         * @param parent
         * @param meshDataContainer
         * @param initialValue
         */
        template<typename _DataType, unsigned int ... _Dimensions>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<_DataType, _Dimensions...>& meshDataContainer,
                                  const DataType& initialValue) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<pos>()>().size(),
                        initialValue);
            Allocator<pos - 1>::allocateMemory(parent, meshDataContainer, initialValue);
        }
    };


    /**
     * @brief The Allocator<_Tp1, dummy> struct
     * Specialization terminating the allocation cycle.
     */
    template<typename dummy>
    struct Allocator<0, dummy>{

        /**
         * @brief Resizes vector on positions according to the respective dimension
         * of mesh and initializes with defaul value if possible.
         * Terminal step for @a pos = 0.
         * @param parent
         * @param mesh
         */
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {

            parent.template getDataByPos<0>().resize(
                        mesh.template getElements<parent.template dimensionAt<0>()>().size());

        }


        /**
         * @brief Resizes vector on positions according to the respective dimension
         * of mesh and initializes with @a initialValue.
         * Terminal step for @a pos = 0.
         * @param parent
         * @param mesh
         * @param initialValue a value to be data initialized with
         */
        template<unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                                  const DataType& initialValue) {

            parent.template getDataByPos<0>().resize(
                        mesh.template getElements<parent.template dimensionAt<0>()>().size(),
                        initialValue);

        }



        /**
         * @brief Allocates memory according to another MashDataContainer which must
         * have data allocated to the @a Dimensions. The data is initialized by
         * default value if possible.
         * @param parent
         * @param meshDataContainer
         */
        template<typename _DataType, unsigned int ... _Dimensions>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<_DataType, _Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<0>()>().size());
        }



        /**
         * @brief Allocates memory according to another %MashDataContainer which must
         * have data allocated to the @a Dimensions. The data is initialized by
         * @a initialValue.
         * @param parent
         * @param meshDataContainer
         * @param initialValue
         */
        template<typename _DataType, unsigned int ... _Dimensions>
        static void allocateMemory(MeshDataContainer<DataType, Dimensions...>& parent ,
                                  const MeshDataContainer<_DataType, _Dimensions...>& meshDataContainer,
                                  const DataType& initialValue) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<0>()>().size(),
                        initialValue);
        }
    };




    /**
     * @brief A structure containing vectors of specified type
     * alocated respectively to the mesh elements.
     */
    _DataContainer<DataType, sizeof... (Dimensions) - 1> data;


public:

    /*
     * Definition of constructors
     */


    MeshDataContainer() = default;

    MeshDataContainer(const MeshDataContainer<DataType, Dimensions...>& meshDataContainer) = default;

    MeshDataContainer(MeshDataContainer<DataType, Dimensions...>&& meshDataContainer) = default;




    /**
     * @brief Resizes the MeshDataContainer according to the mesh.
     */
    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        allocateData(mesh);
    }




    /**
     * @brief Resizes the %MeshDataContainer according to the mesh.
     * Moreover, the @a Dimensions and @a DataType can be deduced from
     * the second and the third arguments.
     */
    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      std::integer_sequence<unsigned int,Dimensions...>, DataType){
        allocateData(mesh);
    }





    /**
     * @brief <HR>
     * This constructor resizes the %MeshDataContainer according to the mesh.
     * And initializes the data with an @a initialValue.
     * @param mesh
     * @param initialValue
     */
    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataType& initialValue){
        allocateData(mesh, initialValue);
    }




    /**
     * @brief <HR>
     * This constructor resizes the %MeshDataContainer according to the mesh.
     * And initializes the data with an @a initialValue.
     * Moreover, the @a Dimensions and @a DataType can be deduced from
     * the second and the third arguments.
     * @param mesh
     * @param initialValue
     */
    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    MeshDataContainer(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                      const DataType& initialValue,
                      std::integer_sequence<unsigned int,Dimensions...>,
                      DataType){
        allocateData(mesh, initialValue);
    }




    /*
     * Data allocation according to mesh or another
     * MeshDataContainer with suitable dimensions
     */

    /**
     * @brief Allocates data according to a @a mesh.
     * @param mesh
     */
    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, mesh);
    }

    /**
     * @brief Allocates data according to a @a mesh and
     * initializes the value by the @a initialValue.
     * @param mesh
     * @param initialValue
     */
    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh, const DataType& initialValue){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, mesh,initialValue);
    }


    /**
     * @brief Allocates data according to another MeshDataContainer.
     * The meshDataContainer mush have @a _Dimensions at least
     * dimensions as @a Dimensions.
     * @param meshDataContainer
     */
    template<typename _DataType, unsigned int ... _Dimensions>
    void allocateData(const MeshDataContainer<_DataType, _Dimensions...>& meshDataContainer){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, meshDataContainer);
    }

    /**
     * @brief Allocates data according to another MeshDataContainer.
     * The meshDataContainer mush have @a _Dimensions at least
     * dimensions as @a Dimensions.
     * @param meshDataContainer
     * @param initialValue
     */
    template<typename _DataType, unsigned int ... _Dimensions>
    void allocateData(const MeshDataContainer<_DataType, _Dimensions...>& meshDataContainer, const DataType& initialValue){
        Allocator<sizeof... (Dimensions) - 1>::allocateMemory(*this, meshDataContainer,initialValue);
    }




    /*
     * Data access
     */

    /**
     * @brief Returns first data array allocated to the given dimension.
     * @return
     */
    template<unsigned int dim>
    DataContainer<DataType, dim>& getDataByDim(){
        return data._DataContainer<DataType, dimensionIndex<dim>()>::_data;
    }

    /**
     * @brief Returns first data array allocated to the given dimension.
     * @return
     */
    template<unsigned int dim>
    const DataContainer<DataType, dim>& getDataByDim() const {
        return data._DataContainer<DataType, dimensionIndex<dim>()>::_data;
    }


    /**
     * @brief getDataByPos Returns data at given position.
     * @return DataContainer at given position.
     */
    template<unsigned int pos>
    DataContainer<DataType, dimensionAt<pos>()>& getDataByPos(){
        return data._DataContainer<DataType,pos>::_data;
    }

    /**
     * @brief getDataByPos
     * @return
     */
    template<unsigned int pos>
    const DataContainer<DataType, dimensionAt<pos>()>& getDataByPos() const {
        return data._DataContainer<DataType,pos>::_data;
    }


    /**
     * @brief Returns an element of the array allocated to the dimension
     * of the passed @a element.
     * @throws out of range
     * @param element A MeshElement of given ElementDimension to select the appropriate data.
     * @return
     */
    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    /**
     * @brief Returns an element of the array allocated to the dimension
     * of the passed @a element.
     * @param element Element of a mesh the data are mapped to
     */
    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    const DataType& at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) const {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    /**
     * @brief Returns the array allocated to position given
     * by the passed integer sequence.
     * @example Example of approaching the arrays by operator[] accepting an integer_sequence
     * @code
     * MeshDataContainer<int, 1,2,3> continer(mesh);
     * constexpr std::integer_sequence< unsigned int, 2 > thirdArray;
     * container[thirdArray][0]; // the first element of the second array
     * @endcode
     */
    template <unsigned int Pos>
    DataContainer<DataType, dimensionAt<Pos>()>& operator[](std::integer_sequence<unsigned int, Pos>) {
        return getDataByPos<Pos>();
    }

    /**
     * @brief Returns the array allocated to position given
     * by the passed integer sequence.
     * @example Example of approaching the arrays by operator[] accepting an integer_sequence
     * @code
     * MeshDataContainer<int, 1,2,3> continer(mesh);
     * constexpr std::integer_sequence< unsigned int, 2 > thirdArray;
     * container[thirdArray][0]; // the first element of the second array
     * @endcode
     */
    template <unsigned int Pos>
    const DataContainer<DataType, dimensionAt<Pos>()>& operator[](std::integer_sequence<unsigned int, Pos>) const {
        return getDataByPos<Pos>();
    }

    /**
     * @brief Returns an element of the array allocated to the dimension
     * of the passed @a element.
     * @param element Element of a mesh the data are mapped to
     */
    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    DataType& operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }

    /**
     * @brief Returns an element of the array allocated to the dimension
     * of the passed @a element.
     * @param element Element of a mesh the data are mapped to
     */
    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    const DataType& operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) const {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }



    /*
     * Assign operators
     */


    /**
     * @brief
     * deep copy operator =
     * @param rhs
     * @return this
     */
    MeshDataContainer<DataType, Dimensions...>& operator=(const MeshDataContainer<DataType, Dimensions...>& rhs) {
        this->data = rhs.data;
        return *this;
    }

    /**
     * @brief
     * move operator =
     * @param rhs
     * @return this
     */
    MeshDataContainer<DataType, Dimensions...>& operator=(MeshDataContainer<DataType, Dimensions...>&& rhs) {
        this->data = std::move(rhs.data);
        return *this;
    }
};





/**
 * @brief The MeshDataContainer<std::tuple<DataTypes>, Dimensions> struct<HR>
 * This specialization of MeshDataContainer allows to declare
 * a data type for each dimension separately. The data types are
 * given by a tuple.
 *
 * @example Definition of MeshDataContainer mapping different types to different dimensions
 * of the mesh.
 * @code
 *
 * MeshDataContainer<std::tuple<int, double>, 3,2> data;
 *
 * @endcode
 */
template <typename ...DataTypes, unsigned int ...Dimensions>
struct MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>{

public:
    template<unsigned int dim>
    static constexpr unsigned int dimensionIndex(){
        return DimensionPos<dim, 0, Dimensions...>::value;
    }

    template<unsigned int pos>
    static constexpr unsigned int dimensionAt(){
        return std::get<pos>(std::array<unsigned int, sizeof... (Dimensions)>{Dimensions...});
    }
public:

    template<unsigned int pos>
    using DataType = typename std::tuple_element<pos,std::tuple<DataTypes...>>::type;
private:
    template<unsigned int Pos, typename Dummy = void>
    struct _DataContainer : _DataContainer<Pos - 1, Dummy>{
        DataContainer<DataType<Pos>, dimensionAt<Pos>()> _data;
    };

    template<typename Dummy>
    struct _DataContainer<0, Dummy>{
        DataContainer<DataType<0>, dimensionAt<0>()> _data;
    };

    /**
     * Data container type according to pos
     */
public:
    template <unsigned int Pos>
    using DataContainerType = DataContainer<DataType<Pos>, dimensionAt<Pos>()>;

    /**
     * @brief size<HR>
     * Returns the number of vectors contained in the MeshDataContainer.
    * @return
    */
   static constexpr unsigned int size() {
       return sizeof... (Dimensions);
   }
private:
    /// Stops the template recursion of allocateMemory
    template <unsigned int pos, typename ContainterType>
    std::enable_if_t<(pos == sizeof...(Dimensions))>
    allocateMemory(const ContainterType& /*mesh*/) {}

    template <unsigned int pos, unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    std::enable_if_t<(pos < sizeof...(Dimensions))>
    allocateMemory(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh) {
        getDataByPos<pos>().resize(mesh.template getElements<dimensionAt<pos>()>().size());
        allocateMemory<pos + 1>(mesh);
    }

    template<unsigned int pos,
             typename _DataType,
             unsigned int Dimension,
             typename IndexType,
             typename Real,
             unsigned int... Reserve,
             typename... _DataTypes>
    std::enable_if_t<(pos < sizeof...(Dimensions))>
    allocateMemory(
        const MeshElements<Dimension, IndexType, Real, Reserve...> &mesh,
        const _DataType &initialValue,
        const _DataTypes &... values)
    {
        getDataByPos<pos>().resize(mesh.template getElements<dimensionAt<pos>()>().size(),
                                   initialValue);
        allocateMemory<pos + 1>(mesh, values...);
    }

    template <unsigned int pos, typename ..._DataTypes, unsigned int ..._Dimensions>
    std::enable_if_t<(pos < sizeof...(Dimensions))>
    allocateMemory(const MeshDataContainer<_DataTypes..., _Dimensions...>& meshDataContainer) {
        getDataByPos<pos>().resize(meshDataContainer.template getDataByDim<dimensionAt<pos>()>().size());
        allocateMemory<pos + 1>(meshDataContainer);
    }

    template<unsigned int pos,
             typename _DataType,
             typename... _DataTypes,
             unsigned int... _Dimensions,
             typename... _ValDataTypes>
    std::enable_if_t<(pos < sizeof...(Dimensions))>
    allocateMemory(
        const MeshDataContainer<std::tuple<_DataTypes...>, _Dimensions...> &meshDataContainer,
        const _DataType &initialValue,
        const _ValDataTypes &... values)
    {
        getDataByPos<pos>()
            .resize(meshDataContainer.template getDataByDim<dimensionAt<pos>()>().size(),
                    initialValue);
        allocateMemory<pos + 1>(meshDataContainer, values...);
    }
/*
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
                meshDataContainer.template getDataByDim<
                        parent.template dimensionAt<pos>()
                    >().size()
                );

            Allocator<pos + 1, _DataTypes...>::allocateMemory(parent, meshDataContainer);
        }

        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer,
                                  const _DataType& initialValue,
                                  const _DataTypes&... values) {
            parent.template getDataByPos<pos>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<pos>()>().size(),
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


        template<unsigned int ... _Dimensions>
        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent ,
                                  const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<0>()>().size());

        }

        template<unsigned int ... _Dimensions>
        static void allocateMemory(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& parent,
                                  const MeshDataContainer<std::tuple<_DataTypes...>, _Dimensions...>& meshDataContainer,
                                  const _DataType& initialValue,
                                  const _DataTypes&...) {
            parent.template getDataByPos<0>().resize(
                        meshDataContainer.template getDataByDim<parent.template dimensionAt<0>()>().size(),
                        initialValue);

        }
    };
*/



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
    DataContainer<std::tuple_element_t<dimensionIndex<dim>(), std::tuple<DataTypes...>>, dim>&
        getDataByDim(){
        return data._DataContainer<dimensionIndex<dim>()>::_data;
    }

    template<unsigned int dim>
    const DataContainer<std::tuple_element_t<dimensionIndex<dim>(), std::tuple<DataTypes...>>, dim>&
        getDataByDim() const {
        return data._DataContainer<dimensionIndex<dim>()>::_data;
    }


    /**
     * @brief GetDataPos
     * @return
     */
    template<unsigned int pos>
    DataContainer<DataType<pos>, dimensionAt<pos>()>& getDataByPos(){
        return data._DataContainer<pos>::_data;
    }

    template<unsigned int pos>
    const DataContainer<DataType<pos>, dimensionAt<pos>()>& getDataByPos() const {
        return data._DataContainer<pos>::_data;
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>&
    at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    const std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>&
    at(const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) const {
        return getDataByDim<ElementDim>().at(element.getIndex());
    }

    /**
     * @brief Returns the array allocated to position given
     * by the passed integer sequence.
     * @example Example of approaching the arrays by operator[] accepting an integer_sequence
     * @code
     * MeshDataContainer<int, 1,2,3> continer(mesh);
     * constexpr std::integer_sequence< unsigned int, 2 > thirdArray;
     * container[thirdArray][0]; // the first element of the second array
     * @endcode
     */
    template <unsigned int Pos>
    DataContainer<std::tuple_element_t<dimensionIndex<Pos>(), std::tuple<DataTypes...>>, Pos>&
    operator[](std::integer_sequence<unsigned int, Pos>) {
        return getDataByPos<Pos>();
    }

    /**
     * @brief Returns the array allocated to position given
     * by the passed integer sequence.
     * @example Example of approaching the arrays by operator[] accepting an integer_sequence
     * @code
     * MeshDataContainer<int, 1,2,3> continer(mesh);
     * constexpr std::integer_sequence< unsigned int, 2 > thirdArray;
     * container[thirdArray][0]; // the first element of the second array
     * @endcode
     */
    template <unsigned int Pos>
    const DataContainer<std::tuple_element_t<dimensionIndex<Pos>(), std::tuple<DataTypes...>>, Pos>&
    operator[](std::integer_sequence<unsigned int, Pos>) const {
        return getDataByPos<Pos>();
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>&
    operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }

    template <unsigned int ElementDim, unsigned int Dimension, typename IndexType, typename Real, unsigned int Reserve>
    const std::tuple_element_t<dimensionIndex<ElementDim>(), std::tuple<DataTypes...>>&
    operator[](const MeshElement<Dimension, ElementDim, IndexType, Real, Reserve>& element) const {
        return getDataByDim<ElementDim>()[element.getIndex()];
    }


    MeshDataContainer() = default;

    MeshDataContainer(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer) = default;

    MeshDataContainer(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>&& meshDataContainer) = default;


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
        allocateMemory<0>(mesh);
    }



    template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
    void allocateData(const MeshElements<Dimension, IndexType, Real, Reserve...>& mesh,
                     const DataTypes&... initialValues){
        allocateMemory<0>(mesh, initialValues...);
    }


    void allocateData(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer){
        allocateMemory<0>(meshDataContainer);
    }



    void allocateData(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& meshDataContainer,
                      const DataTypes&... initialValues){
        allocateMemory<0>(meshDataContainer, initialValues...);
    }

    MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& operator=(const MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& rhs) {
        this->data = rhs.data;
        return *this;
    }

    MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>& operator=(MeshDataContainer<std::tuple<DataTypes...>, Dimensions...>&& rhs) {
        this->data = std::move(rhs.data);
        return *this;
    }

};




/**
 * @brief MakeMeshDataContainer<HR>
 * struct publishing the type of %MeshDataContainer using
 * tuple and integer sequence
 */
template<typename DataType, typename Sequence>
struct MakeMeshDataContainer {
    static_assert(
        std::is_class<Sequence>::value && !std::is_class<Sequence>::value,
        "The Sequence parameter in MakeMeshDataContainer "
        "must be a std::integer_sequence<unsigned int, seq...>"
        "please notice the type unsigned int"
        );
};

template<typename Type, unsigned int... Dimensions>
struct MakeMeshDataContainer<Type, std::integer_sequence<unsigned int, Dimensions...>>{
    using type = MeshDataContainer<Type, Dimensions...>;
};

template<typename DataType, typename Sequence>
using MakeMeshDataContainer_t = typename MakeMeshDataContainer<DataType, Sequence>::type;

#endif // MESHDATACONTAINER_H
