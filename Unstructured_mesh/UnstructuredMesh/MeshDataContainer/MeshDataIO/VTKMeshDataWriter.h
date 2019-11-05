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
    static void writeColumn(...){}

    template<typename T, unsigned int Index, unsigned int Position, typename IndexType, typename Real>
    static auto writeColumn(std::ostream& ost, const DataContainer<T, Position, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer)
    -> typename std::enable_if<
        Detail::is_indexable<typename Traits<T>::ttype::template type<Index>>::value
    >::type
    {
        if (data.at(0).size() == MeshDimension)
        ost << "VECTORS ";
        else if (data.at(0).size() == MeshDimension * MeshDimension)
        ost << "TENZORS ";

        ost << Traits<T>::ttype::template getName<Index>() << "double";

        IndexType i = 0;
        for (; i < data.size(); i++) {
            for (unsigned int j = 0; j < data.at(i).size(); j++) {
                ost << data.at(i)[j];
            }

        }
    }



    template<typename T,unsigned int Index = 0, typename VOID = void>
    struct writeCellData{};

    template<typename T,unsigned int Index, typename... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index < Traits<T, Types...>::size() - 1>>{
        template<unsigned int Position, typename IndexType, typename Real>
        static void write(std::ostream& ost, const DataContainer<T, Position, MeshDimension> &data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer){

            ost << "CELL_DATA " << writer.cellVert.template getDataByPos<0>().size() << std::endl;

            ost << '"' << Traits<T, Types...>::template getName<Index>() << "\" : ";
            VariableExport::_writeWar(ost, Traits<T, Types...>::template getReference<Index>()->getValue(data));
            ost << ", ";
            writeCellData<Traits<T, Types...>, Index + 1>::print(ost, data);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct writeCellData <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{
        static void write(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T, Types...>::template getName<Traits<T, Types...>::size() - 1>() << "\" : ";
            VariableExport::_writeWar(ost, Traits<T, Types...>::template getReference<Traits<T, Types...>::size() - 1>()->getValue(traitedClass));
        }
    };




public:
    template<typename T,typename IndexType, typename Real, unsigned int ...Dimensions>
    static void writeToStream(std::ostream& ost, MeshDataContainer<T, Dimensions...>& data, VTKMeshWriter<MeshDimension,IndexType, Real>& writer) {
        using type = typename decltype(data.template getDataByDim<MeshDimension>())::DataType;
        static_assert (Detail::has_default_traits<type>::value, "The class T must have defined traits for example using macro MAKE_ATTRIBUTE_TRAIT in header Traits.h");

    }

};


#endif // VTKMESHDATAWRITER_H
