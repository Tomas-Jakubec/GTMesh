#ifndef VTKMESHDATAWRITER_H
#define VTKMESHDATAWRITER_H
#include "../../../Traits/Traits.h"
#include "../../../Traits/CustomTypeTraits.h"
#include "../MeshDataContainer.h"
#include "../../../Debug/Debug.h"
#include "../../MeshIO/MeshWriter/VTKMeshWriter.h"

#include <ostream>

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
    static auto writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    -> typename std::enable_if<
        IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
        MeshDimension == 3
       >::type
    {

        if (DefaultIOTraits<T>::tr.template getValue<Index>(data.at(0)).size() == MeshDimension)
            ost << "VECTORS ";
        else if (DefaultIOTraits<T>::tr.template getValue<Index>(data.at(0)).size() == MeshDimension * MeshDimension)
            ost << "TENZORS ";

        ost << DefaultIOTraits<T>::tr.template getName<Index>() << " double\n";


        IndexType realIndex = 0;
        IndexType localIndex = 0;
        for(const std::pair<IndexType, IndexType>& key : writer.backwardCellIndexMapping) {
            while (localIndex < key.first) {
                    for (unsigned int j = 0; j < DefaultIOTraits<T>::tr.template getValue<Index>(data.at(0)).size(); j++) {
                    ost << DefaultIOTraits<T>::tr.template getValue<Index>(data.at(realIndex))[j] << ' ';
                }
                realIndex++;
                localIndex++;
            }
            realIndex = key.second;
            localIndex++;
            for (unsigned int j = 0; j < DefaultIOTraits<T>::tr.template getValue<Index>(data.at(0)).size(); j++) {
                ost << DefaultIOTraits<T>::tr.template getValue<Index>(data.at(realIndex))[j] << ' ';
            }
        }
        while (realIndex < data.size() - 1) {
            for (unsigned int j = 0; j < DefaultIOTraits<T>::tr.template getValue<Index>(data.at(0)).size(); j++) {
                ost << DefaultIOTraits<T>::tr.template getValue<Index>(data.at(realIndex))[j] << ' ';
            }
            realIndex++;
        }
    }




    template<typename T, unsigned int Index, typename IndexType, typename Real>
    static auto writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    -> typename std::enable_if<
        IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
        MeshDimension == 2
       >::type
    {

        if (DefaultIOTraits<T>::traitsType::template getReference<Index>()->getValue(data.at(0)).size() == MeshDimension)
            ost << "VECTORS ";
        else if (DefaultIOTraits<T>::traitsType::template getReference<Index>()->getValue(data.at(0)).size() == MeshDimension * MeshDimension)
            ost << "TENZORS ";

        ost << DefaultIOTraits<T>::traitsType::template getName<Index>() << " double\n";


        IndexType realIndex = 0;
        IndexType localIndex = 0;
        for(const std::pair<IndexType, IndexType>& key : writer.backwardCellIndexMapping) {
            while (localIndex < key.first) {
                    for (unsigned int j = 0; j < DefaultIOTraits<T>::traitsType::template getReference<Index>()->getValue(data.at(0)).size(); j++) {
                    ost << DefaultIOTraits<T>::traitsType::template getValue<Index>(data.at(realIndex))[j] << " 0.0 ";
                }
                realIndex++;
                localIndex++;
            }
            realIndex = key.second;
            localIndex++;
            for (unsigned int j = 0; j < DefaultIOTraits<T>::traitsType::template getReference<Index>()->getValue(data.at(0)).size(); j++) {
                ost << DefaultIOTraits<T>::traitsType::template getValue<Index>(data.at(realIndex))[j] << " 0.0 ";
            }
        }
        while (realIndex < data.size() - 1) {
            for (unsigned int j = 0; j < DefaultIOTraits<T>::traitsType::template getReference<Index>()->getValue(data.at(0)).size(); j++) {
                ost << DefaultIOTraits<T>::traitsType::template getValue<Index>(data.at(realIndex))[j] << " 0.0 ";
            }
            realIndex++;
        }
    }

    template<typename T, unsigned int Index, typename IndexType, typename Real>
    static auto writeColumn(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    -> typename std::enable_if<
        !IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value
    >::type
    {


        ost << "SCALARS " << DefaultIOTraits<T>::tr.template getName<Index>() << " double 1\nLOOKUP_TABLE default\n";

        IndexType realIndex = 0;
        IndexType localIndex = 0;
        for(const std::pair<IndexType, IndexType>& key : writer.backwardCellIndexMapping) {
            while (localIndex < key.first) {
                ost << DefaultIOTraits<T>::tr.template getValue<Index>(data.at(realIndex)) << ' ';
                realIndex++;
                localIndex++;
            }
            realIndex = key.second;
            localIndex++;
            ost << DefaultIOTraits<T>::tr.template getValue<Index>(data.at(realIndex)) << ' ';
        }
        while (realIndex < data.size() - 1) {
            ost << DefaultIOTraits<T>::tr.template getValue<Index>(data.at(realIndex)) << ' ';
            realIndex++;
        }
    }


    template<typename T,unsigned int Index = 0, typename Void = void>
    struct writeCellData{};

    template<typename T,unsigned int Index, typename... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index < Traits<T, Types...>::size() - 1>>{

        template<typename IndexType, typename Real>

        static void write(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){
            DBGVAR(IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value);
            writeColumn<T, Index, IndexType, Real>(ost, data, writer);
            ost << std::endl;
            writeCellData<Traits<T, Types...>, Index + 1>::write(ost, data, writer);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{
        template< typename IndexType, typename Real>
        static void write(std::ostream& ost, const DataContainer<T, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){

            writeColumn<T, Index, IndexType, Real>(ost, data, writer);
            ost << std::endl;
        }
    };


public:
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
        static auto writeToStream(std::ostream& ost,
                                  MeshDataContainer<T, Dimensions...>& data,
                                  VTKMeshWriter<MeshDimension,IndexType, Real>& writer) -> typename
        std::enable_if<
            (!HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type>::value) ||
            (!(MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::getMappedDimension() == MeshDimension))
        >::type
        {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type;


            MeshDataIterator<Index - 1, OK | HasDefaultIOTraits<type>::value>:: writeToStream(ost, data, writer);
        }

        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static auto writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
            ->typename std::enable_if<HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type>::value &&
                              MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::getMappedDimension() == MeshDimension
                             >::type
        {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type;

            VTKMeshDataWriter<MeshDimension>::writeCellData<typename DefaultIOTraits<type>::traitsType>::write(ost, data.template getDataByPos<Index>(), writer);

            MeshDataIterator<Index - 1, OK | HasDefaultIOTraits<type>::value>:: writeToStream(ost, data, writer);
        }
    };

    template <bool OK>
    struct MeshDataIterator <0, OK> {
        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static auto writeToStream(std::ostream&,
                                  MeshDataContainer<T, Dimensions...>&,
                                  VTKMeshWriter<MeshDimension,IndexType, Real>&) -> typename
            std::enable_if<
                (!HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type>::value) ||
                (!(MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::getMappedDimension() == MeshDimension))
            >::type
        {
            static_assert (OK , "The mesh data container must have at least one DataContainer mapped to cells with traits for example using macro MAKE_ATTRIBUTE_TRAIT see header Traits.h");

        }

        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static auto writeToStream(std::ostream& ost,
                                  MeshDataContainer<T, Dimensions...>& data,
                                  VTKMeshWriter<MeshDimension,IndexType, Real>& writer) -> typename
            std::enable_if<
                HasDefaultIOTraits<typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type>::value &&
                MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::getMappedDimension() == MeshDimension
            >::type
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
