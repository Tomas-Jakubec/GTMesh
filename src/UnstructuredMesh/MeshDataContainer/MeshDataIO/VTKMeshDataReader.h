#ifndef VTKMESHDATAREADER_H
#define VTKMESHDATAREADER_H
#include "../../../Traits/Traits.h"
#include "../../../Traits/CustomTypeTraits.h"
#include "../MeshDataContainer.h"
#include "../../../Debug/Debug.h"
#include <istream>
#include <map>
#include <sstream>

template <unsigned int MeshDimension, typename IndexType>
class VTKMeshDataReader {

    static_assert (MeshDimension == 2 || MeshDimension == 3, "The VTK file format can represent data only in 2D or 3D");
    /**
     * @brief readColumn
     * reads a single column of traited data
     */
    static void readColumn(std::istream& ist [[maybe_unused]],...){
        DBGMSG("capture");
        throw std::runtime_error("capture of read column must not be called.");
    }

    template<typename T, unsigned int Index, unsigned int Position>
    static auto readColumn(std::istream& ist, DataContainer<T, Position, MeshDimension> &data,std::map<std::string, std::istream::pos_type>& dataPositions)
    -> typename std::enable_if<
        IsIndexable<typename Traits<T>::ttype::template type<Index>>::value &&
        MeshDimension == 3
       >::type
    {

        ist.seekg(dataPositions[Traits<T>::ttype::template getName<Index>()]);

        typename Traits<T>::ttype::template type<Index> value;

        for (IndexType i = 0; i < data.size(); i++) {
            for (unsigned int j = 0; j < Traits<T>::ttype::template getValue<Index>(data.at(i)).size(); j++){
                ist >> value[j];
            }
            Traits<T>::ttype::template setValue<Index>(data.at(i), value);
        }

    }


    template<typename T, unsigned int Index, unsigned int Position>
    static auto readColumn(std::istream& ist, DataContainer<T, Position, MeshDimension> &data,std::map<std::string, std::istream::pos_type>& dataPositions)
    -> typename std::enable_if<
        IsIndexable<typename Traits<T>::ttype::template type<Index>>::value &&
        MeshDimension == 2
       >::type
    {

        ist.seekg(dataPositions[Traits<T>::ttype::template getName<Index>()]);

        typename Traits<T>::ttype::template type<Index> value;
        typename Traits<T>::ttype::template type<Index> dummy;

        for (IndexType i = 0; i < data.size(); i++) {
            for (unsigned int j = 0; j < Traits<T>::ttype::template getValue<Index>(data.at(i)).size(); j++){
                ist >> value[j];
            }

            ist >> dummy[0];

            Traits<T>::ttype::template setValue<Index>(data.at(i), value);
        }

    }


    template<typename T, unsigned int Index, unsigned int Position>
    static auto readColumn(std::istream& ist, DataContainer<T, Position, MeshDimension> &data,std::map<std::string, std::istream::pos_type>& dataPositions)
    -> typename std::enable_if<
        !IsIndexable<typename Traits<T>::ttype::template type<Index>>::value
    >::type
    {


        ist.seekg(dataPositions[Traits<T>::ttype::template getName<Index>()]);

        typename Traits<T>::ttype::template type<Index> value;

        for (IndexType i = 0; i < data.size(); i++){
            ist >> value;
            Traits<T>::ttype::template setValue<Index>(data.at(i), value);
        }

    }
private:

    template<typename T,unsigned int Index = 0, typename Void = void>
    struct readCellData{};

    template<typename T,unsigned int Index, typename... Types>
    struct readCellData <Traits<T, Types...>, Index, std::enable_if_t<Index < Traits<T, Types...>::size() - 1>>{

        template<unsigned int Position>

        static void read(std::istream& ist, DataContainer<T, Position, MeshDimension> &data, std::map<std::string, std::istream::pos_type>& dataPositions){
            DBGVAR(IsIndexable<typename Traits<T>::ttype::template type<Index>>::value);
            readColumn<T, Index, Position>(ist, data, dataPositions);
            readCellData<Traits<T, Types...>, Index + 1>::read(ist, data, dataPositions);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct readCellData <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{
        template<unsigned int Position>
        static void read(std::istream& ist, DataContainer<T, Position, MeshDimension> &data, std::map<std::string, std::istream::pos_type>& dataPositions){

            readColumn<T, Index, Position>(ist, data, dataPositions);

        }
    };



public:

    static std::map<std::string, std::istream::pos_type> indexData(std::istream& ist) {

        std::map<std::string, std::istream::pos_type> dataPositions;

        std::string line;
        ist.seekg(ist.beg);
        while(getline(ist, line)) {
            int flag = (line.find("SCALARS")!= line.npos ? 1
                             : line.find("VECTORS") != line.npos ? 2
                                   : 0 );
            if (flag != 0){
                std::string dataName;
                std::stringstream sstream(line);
                sstream.ignore(9, ' ');
                sstream >> dataName;

                if (flag == 1) { // scalar quantity found
                    ist.ignore(500, '\n');
                }

                dataPositions.insert(std::pair(dataName, ist.tellg()));
            }
        }
        ist.clear();
        return dataPositions;
    }

    template<typename T, unsigned int Position>
    static void readData(std::istream& ist, DataContainer<T, Position, MeshDimension>& data) {

        std::map<std::string, std::istream::pos_type> dataPositions = indexData(ist);

        readCellData<typename Traits<T>::ttype>::read(ist, data, dataPositions);
    }
};


#endif // VTKMESHDATAREADER_H
