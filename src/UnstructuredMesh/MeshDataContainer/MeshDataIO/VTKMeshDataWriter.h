#ifndef VTKMESHDATAWRITER_H
#define VTKMESHDATAWRITER_H
#include "../../../Traits/Traits.h"
#include "../../../Traits/CustomTypeTraits.h"
#include "../MeshDataContainer.h"
#include "../../../Debug/Debug.h"
#include "../../MeshIO/MeshWriter/VTKMeshWriter.h"

#include <ostream>
/**
 * @brief The VTKMeshDataWriter utilizes the
 * class Traits to export the given MeshDataContainer
 * of DataContainer automatically. If the MeshDataContainer
 * is given, then it is iterated and the function looks
 * for the exportable data structures with DefaultIOTraits.
 * @example
 * @code
 * struct DataStruct{
 *    double pressure;
 *    Vector<double, 3> velocity;
 * }
 * MAKE_ATTRIBUTE_TRAIT(DataStruct, pressure, velocity);
 * // Allocate the data with pressure = 42, velocity = {1,2,3}
 * MeshDataContainer<DataStruct,3> meshData(mesh, {42, {1,2,3}});
 * // Mesh export etc.
 * VTKMeshDataWriter::wrtieToStream(cerr, meshData, meshWriter);
 * // Writes "CELL_DATA %number_of_cells
 * // SCALAR presure
 * // 42 42 ...
 * // VECTOR velocity
 * // 1 2 3 1 2 3 ..."
 * @endcode
 */
template <unsigned int MeshDimension>
class VTKMeshDataWriter {

    static_assert (MeshDimension == 2 || MeshDimension ==3, "The VTK format can hold only 2D or 3D data");

    /**
     * @brief writeColumn
     * writes a single column of traited data
     */
    static void writeColumn(std::ostream&,...){
        DBGMSG("capture");
        throw std::runtime_error("capture of write column must not be called.");
    }

    template<typename T, unsigned int Index, typename IndexType, typename Real>
    static
    typename std::enable_if<
       IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
       MeshDimension == 3
    >::type
    writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    {


        if (DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size() == MeshDimension)
            ost << "VECTORS ";
        else if (DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size() == MeshDimension * MeshDimension)
            ost << "TENZORS ";

        ost << DefaultIOTraits<T>::getTraits().template getName<Index>() << " double\n";


        IndexType realIndex = 0;
        IndexType localIndex = 0;
        auto keyIt = writer.backwardCellIndexMapping.cbegin();
        while(keyIt != writer.backwardCellIndexMapping.cend()) {
            const std::pair<IndexType, IndexType>& key = *keyIt;
            while (localIndex < key.first) {
                    for (unsigned int j = 0; j < DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size(); j++) {
                    ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex))[j] << ' ';
                }
                realIndex++;
                localIndex++;
            }
            realIndex = key.second;
            localIndex++;
            for (unsigned int j = 0; j < DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size(); j++) {
                ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex))[j] << ' ';
            }
            ++keyIt;
            if (keyIt == writer.backwardCellIndexMapping.cend() || realIndex != keyIt->second){
                realIndex++;
            }
        }

        // write the rest of not tessellated cells
        while (realIndex < data.size()) {
            for (unsigned int j = 0; j < DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size(); j++) {
                ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex))[j] << ' ';
            }
            realIndex++;
        }
    }

    template<typename T, unsigned int Index, typename IndexType, typename Real>
    static
    typename std::enable_if<
       !IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
       MeshDimension == 3
    >::type
    writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    {


        ost << "SCALARS " << DefaultIOTraits<T>::getTraits().template getName<Index>() << " double 1\nLOOKUP_TABLE default\n";

        IndexType realIndex = 0;
        IndexType localIndex = 0;
        auto keyIt = writer.backwardCellIndexMapping.cbegin();
        while(keyIt != writer.backwardCellIndexMapping.cend()) {
            const std::pair<IndexType, IndexType>& key = *keyIt;
            while (localIndex < key.first) {
                ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex)) << ' ';
                realIndex++;
                localIndex++;
            }
            realIndex = key.second;
            localIndex++;
            ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex)) << ' ';
            ++keyIt;
            if (keyIt == writer.backwardCellIndexMapping.cend() || realIndex != keyIt->second){
                realIndex++;
            }
        }
        while (realIndex < data.size()) {
            ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex)) << ' ';
            realIndex++;
        }
    }




    template<typename T, unsigned int Index, typename IndexType, typename Real>
    static
    typename std::enable_if<
       IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
       MeshDimension == 2
    >::type
    writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>&)
    {

        if (DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size() == MeshDimension)
            ost << "VECTORS ";
        else if (DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size() == MeshDimension * MeshDimension)
            ost << "TENZORS ";

        ost << DefaultIOTraits<T>::getTraits().template getName<Index>() << " double\n";

        for (const auto& d : data){
            for(unsigned int j = 0; j < DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(0)).size(); j++){
                ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(d)[j] << ' ';
            }
            ost << "0.0 ";
        }
    }

    template<typename T, unsigned int Index, typename IndexType, typename Real>
    static
    typename std::enable_if<
        !IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
        MeshDimension == 2
    >::type
    writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>&)
    {

        ost << "SCALARS " << DefaultIOTraits<T>::getTraits().template getName<Index>() << " double 1\nLOOKUP_TABLE default\n";

        IndexType realIndex = 0;
        while (realIndex < data.size()) {
            ost << DefaultIOTraits<T>::getTraits().template getValue<Index>(data.at(realIndex)) << ' ';
            realIndex++;
        }
    }




    template<typename T,unsigned int Index = 0, typename Void = void>
    struct writeCellData{};

    template<typename T,unsigned int Index, typename... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<(Index < Traits<T, Types...>::size() - 1)>>{

        template<typename IndexType, typename Real>

        static void write(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){
            //DBGVAR(IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value);
            writeColumn<T, Index, IndexType, Real>(ost, data, writer);
            ost << std::endl;
            writeCellData<Traits<T, Types...>, Index + 1>::write(ost, data, writer);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{
        template< typename IndexType, typename Real>
        static void write(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){
            //DBGVAR(IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value);
            writeColumn<T, Index, IndexType, Real>(ost, data, writer);
            ost << std::endl;
        }
    };


public:
    /**
     * @brief Exports the given data in the given DataContainer in VTK format
     * to the ostream ost. The given data type T must have defined DefaultIOTraits.
     * To define the traits, one can utilize the macro MAKE_ATTRIBUTE_TRAIT(_IO), @see Traits.
     * The data are axported with the line CELL_DATA.
     * @example
     * @code
     * struct DataStruct{
     *    double pressure;
     *    Vector<double, 3> velocity;
     * }
     * MAKE_ATTRIBUTE_TRAIT(DataStruct, pressure, velocity);
     * // Allocate the data with pressure = 42, velocity = {1,2,3}
     * MeshDataContainer<DataStruct,3> meshData(mesh, {42, {1,2,3}});
     * // Mesh export etc.
     * VTKMeshDataWriter::wrtieToStream(cerr, meshData.getDataByPos<0>(), meshWriter);
     * @endcode
     */
    template<typename T,typename IndexType, typename Real>
    static void writeToStream(std::ostream& ost, DataContainer<T, MeshDimension>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
        using type = T;//typename std::remove_reference<decltype(data.template getDataByDim<MeshDimension>())>::type::DataType;
        static_assert (HasDefaultIOTraits<type>::value, "The class T must have defined traits for example using macro MAKE_ATTRIBUTE_TRAIT in header Traits.h");
        ost << "CELL_DATA " << writer.getNumberOfCells() << std::endl;
        writeCellData<typename DefaultIOTraits<type>::traitsType>::write(ost, data, writer);
    }

private:


    template <unsigned int Index, bool OK = false>
    struct MeshDataIterator{

        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static
        typename std::enable_if<
            (!HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type>::value) ||
            (!(MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::getMappedDimension() == MeshDimension))
        >::type
        writeToStream(std::ostream& ost,
                                  MeshDataContainer<T, Dimensions...>& data,
                                  VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
        {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type;


            MeshDataIterator<Index - 1, OK | HasDefaultIOTraits<type>::value>:: writeToStream(ost, data, writer);
        }

        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static
        typename std::enable_if<
            HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type>::value &&
            MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::getMappedDimension() == MeshDimension
        >::type
        writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)

        {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type;

            VTKMeshDataWriter<MeshDimension>::writeCellData<typename DefaultIOTraits<type>::traitsType>::write(ost, data.template getDataByPos<Index>(), writer);

            MeshDataIterator<Index - 1, OK | HasDefaultIOTraits<type>::value>:: writeToStream(ost, data, writer);
        }
    };

    template <bool OK>
    struct MeshDataIterator <0, OK> {
        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static
        typename std::enable_if<
            (!HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type>::value) ||
            (!(MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::getMappedDimension() == MeshDimension))
        >::type
        writeToStream(std::ostream&,
                                  MeshDataContainer<T, Dimensions...>&,
                                  VTKMeshWriter<MeshDimension,IndexType, Real>&)
        {
            static_assert (OK , "The mesh data container must have at least one DataContainer mapped to cells with traits for example using macro MAKE_ATTRIBUTE_TRAIT see header Traits.h");

        }

        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static
        typename std::enable_if<
            HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type>::value &&
            MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::getMappedDimension() == MeshDimension
        >::type
        writeToStream(std::ostream& ost,
                                  MeshDataContainer<T, Dimensions...>& data,
                                  VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
        {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type;

            VTKMeshDataWriter<MeshDimension>::writeCellData<typename DefaultIOTraits<type>::traitsType>::write(ost, data.template getDataByPos<0>(), writer);

        }
    };

public:
    template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
    static void writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
        using type = T;//typename std::remove_reference<decltype(data.template getDataByDim<MeshDimension>())>::type::DataType;
        static_assert (HasDefaultIOTraits<type>::value, "The class T must have defined traits for example using macro MAKE_ATTRIBUTE_TRAIT in header Traits.h");
        ost << "CELL_DATA " << writer.getNumberOfCells() << std::endl;
        MeshDataIterator<MeshDataContainer<T, Dimensions...>::size() - 1>::writeToStream(ost, data, writer);
    }
};

#endif // VTKMESHDATAWRITER_H
