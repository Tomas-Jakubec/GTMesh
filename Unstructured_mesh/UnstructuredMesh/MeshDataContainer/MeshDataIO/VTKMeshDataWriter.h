#ifndef VTKMESHDATAWRITER_H
#define VTKMESHDATAWRITER_H
#include "../Traits.h"
#include "../MeshDataContainer.h"
#include "../../MeshIO/MeshWriter/VTKMeshWriter.h"

#include <ostream>

template <unsigned int MeshDimension>
class VTKMeshDataWriter {


    /**
     * @brief writeColumn
     * writes a single column of traited data
     */
    static void writeColumn(std::ostream& ost [[maybe_unused]],...){
        DBGMSG("capture");
        throw std::runtime_error("capture of write column must not be called.");
    }

    template<typename T, unsigned int Index, unsigned int Position, typename IndexType, typename Real>
    static auto writeColumn(std::ostream& ost, const DataContainer<T, Position, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    -> typename std::enable_if<
        Detail::is_indexable<typename Traits<T>::ttype::template type<Index>>::value
       >::type
    {

        if (Traits<T>::ttype::template getReference<Index>()->getValue(data.at(0)).size() == MeshDimension)
        ost << "VECTORS ";
        else if (Traits<T>::ttype::template getReference<Index>()->getValue(data.at(0)).size() == MeshDimension * MeshDimension)
        ost << "TENZORS ";

        ost << Traits<T>::ttype::template getName<Index>() << " double\n";


        IndexType realIndex = 0;
        for (IndexType i = 0; i < writer.getNumberOfCells(); i++) {

            auto iterator = writer.backwardCellIndexMapping.find(i);
            if (iterator == writer.backwardCellIndexMapping.end()){
                for (unsigned int j = 0; j < Traits<T>::ttype::template getReference<Index>()->getValue(data.at(0)).size(); j++) {
                    ost << Traits<T>::ttype::template getValue<Index>(data.at(realIndex))[j] << ' ';
                }
                realIndex++;
            } else {
                realIndex = iterator->second;
                for (unsigned int j = 0; j < Traits<T>::ttype::template getReference<Index>()->getValue(data.at(0)).size(); j++) {
                    ost << Traits<T>::ttype::template getValue<Index>(data.at(realIndex))[j] << ' ';
                }
            }
        }
    }


    template<typename T, unsigned int Index, unsigned int Position, typename IndexType, typename Real>
    static auto writeColumn(std::ostream& ost, const DataContainer<T, Position, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    -> typename std::enable_if<
        !Detail::is_indexable<typename Traits<T>::ttype::template type<Index>>::value
    >::type
    {


        ost << "SCALARS " << Traits<T>::ttype::template getName<Index>() << " double 1\nLOOKUP_TABLE default\n";

        IndexType realIndex = 0;
        for (IndexType i = 0; i < writer.getNumberOfCells(); i++) {

            auto iterator = writer.backwardCellIndexMapping.find(i);
            if (iterator == writer.backwardCellIndexMapping.end()){

                ost << Traits<T>::ttype::template getValue<Index>(data.at(realIndex)) << ' ';

                realIndex++;
            } else {

                realIndex = iterator->second;
                ost << Traits<T>::ttype::template getValue<Index>(data.at(realIndex)) << ' ';

            }
        }
    }


    template<typename T,unsigned int Index = 0, typename VOID = void>
    struct writeCellData{};

    template<typename T,unsigned int Index, typename... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index < Traits<T, Types...>::size() - 1>>{

        template<unsigned int Position, typename IndexType, typename Real>

        static void write(std::ostream& ost, const DataContainer<T, Position, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){
            DBGVAR(Detail::is_indexable<typename Traits<T>::ttype::template type<Index>>::value);
            writeColumn<T, Index, Position, IndexType, Real>(ost, data, writer);
            ost << std::endl;
            writeCellData<Traits<T, Types...>, Index + 1>::write(ost, data, writer);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{
        template<unsigned int Position, typename IndexType, typename Real>
        static void write(std::ostream& ost, const DataContainer<T, Position, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){

            writeColumn<T, Index, Position, IndexType, Real>(ost, data, writer);
            ost << std::endl;
        }
    };


public:
    template<typename T,typename IndexType, typename Real, unsigned int Position>
    static void writeToStream(std::ostream& ost, DataContainer<T, Position, MeshDimension>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
        using type = T;//typename std::remove_reference<decltype(data.template getDataByDim<MeshDimension>())>::type::DataType;
        static_assert (Detail::has_default_traits<type>::value, "The class T must have defined traits for example using macro MAKE_ATTRIBUTE_TRAIT in header Traits.h");
        ost << "CELL_DATA " << writer.getNumberOfCells() << std::endl;
        writeCellData<typename Traits<type>::ttype>::write(ost, data, writer);
    }

private:
    template <unsigned int Index, bool OK = false>
    struct MeshDataIterator{
        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static void writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type;

            if constexpr (Detail::has_default_traits<type>::value && MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::getMappedDimension() == MeshDimension){
                VTKMeshDataWriter<MeshDimension>::writeCellData<typename Traits<type>::ttype>::write(ost, data.template getDataByPos<Index>(), writer);
            }
            MeshDataIterator<Index - 1, OK | Detail::has_default_traits<type>::value>:: writeToStream(ost, data, writer);
        }
    };

    template <bool OK>
    struct MeshDataIterator <0, OK> {
        template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
        static void writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
            using type = typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type;
            static_assert (OK | Detail::has_default_traits<type>::value, "The mesh data container must have at least one DataContainer mapped to cells with traits for example using macro MAKE_ATTRIBUTE_TRAIT see header Traits.h");
            if constexpr (Detail::has_default_traits<type>::value && MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::getMappedDimension() == MeshDimension){
                VTKMeshDataWriter<MeshDimension>::writeCellData<typename Traits<type>::ttype>::write(ost, data.template getDataByPos<0>(), writer);
            }
        }
    };

public:
    template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
    static void writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
        //using type = T;//typename std::remove_reference<decltype(data.template getDataByDim<MeshDimension>())>::type::DataType;
        //static_assert (Detail::has_default_traits<type>::value, "The class T must have defined traits for example using macro MAKE_ATTRIBUTE_TRAIT in header Traits.h");
        ost << "CELL_DATA " << writer.getNumberOfCells() << std::endl;
        MeshDataIterator<MeshDataContainer<T, Dimensions...>::size() - 1>::writeToStream(ost, data, writer);
    }
};


#endif // VTKMESHDATAWRITER_H
