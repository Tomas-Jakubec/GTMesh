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
    static void readColumn(std::istream& ,...){
        DBGMSG("capture");
        throw std::runtime_error("capture of read column must not be called.");
    }

    template<typename T, unsigned int Index>
    static auto readColumn(std::istream& ist, DataContainer<T, MeshDimension> &data,std::map<std::string, std::istream::pos_type>& dataPositions)
    -> typename std::enable_if<
        IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
        MeshDimension == 3
       >::type
    {

        ist.seekg(dataPositions[DefaultIOTraits<T>::tr.template getName<Index>()]);
    std::string line;
    std::getline(ist, line);
        ist.seekg(dataPositions[DefaultIOTraits<T>::tr.template getName<Index>()]);

        typename DefaultIOTraits<T>::traitsType::template type<Index> value;

        for (IndexType i = 0; i < data.size(); i++) {
            for (unsigned int j = 0; j < DefaultIOTraits<T>::tr.template getValue<Index>(data.at(i)).size(); j++){
                ist >> value[j];
            }
            DefaultIOTraits<T>::tr.template setValue<Index>(data.at(i), value);
        }

    }


    template<typename T, unsigned int Index>
    static auto readColumn(std::istream& ist, DataContainer<T, MeshDimension> &data,std::map<std::string, std::istream::pos_type>& dataPositions)
    -> typename std::enable_if<
        IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
        MeshDimension == 2
       >::type
    {

        ist.seekg(dataPositions[DefaultIOTraits<T>::traitsType::template getName<Index>()]);

        typename DefaultIOTraits<T>::traitsType::template type<Index> value;
        typename DefaultIOTraits<T>::traitsType::template type<Index> dummy;

        for (IndexType i = 0; i < data.size(); i++) {
            for (unsigned int j = 0; j < DefaultIOTraits<T>::tr.template getValue<Index>(data.at(i)).size(); j++){
                ist >> value[j];
            }

            ist >> dummy[0];

            DefaultIOTraits<T>::tr.template setValue<Index>(data.at(i), value);
        }

    }


    template<typename T, unsigned int Index>
    static auto readColumn(std::istream& ist, DataContainer<T, MeshDimension> &data,std::map<std::string, std::istream::pos_type>& dataPositions)
    -> typename std::enable_if<
        !IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value
    >::type
    {


        ist.seekg(dataPositions[DefaultIOTraits<T>::tr.template getName<Index>()]);

        typename DefaultIOTraits<T>::traitsType::template type<Index> value;

        for (IndexType i = 0; i < data.size(); i++){
            ist >> value;
            DefaultIOTraits<T>::tr.template setValue<Index>(data.at(i), value);
        }

    }
private:

    template<typename T,unsigned int Index = 0, typename Void = void>
    struct readCellData{};

    template<typename T,unsigned int Index, typename... Types>
    struct readCellData <Traits<T, Types...>, Index, std::enable_if_t<Index < Traits<T, Types...>::size() - 1>>{



        static void read(std::istream& ist, DataContainer<T, MeshDimension> &data, std::map<std::string, std::istream::pos_type>& dataPositions){

            readColumn<T, Index>(ist, data, dataPositions);
            readCellData<Traits<T, Types...>, Index + 1>::read(ist, data, dataPositions);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct readCellData <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{

        static void read(std::istream& ist, DataContainer<T, MeshDimension> &data, std::map<std::string, std::istream::pos_type>& dataPositions){

            readColumn<T, Index>(ist, data, dataPositions);

        }
    };



public:

    static std::map<std::string, std::istream::pos_type> indexData(std::istream& ist) {

        std::map<std::string, std::istream::pos_type> dataPositions;
        /*
        if ((ist.flags() & std::ios::binary) != 0) {
            std::runtime_error("open the file stream as binary to ensure correct behaviour of tellg");
        }
        */

        std::string line;
        ist.seekg(ist.beg);
        while(std::getline(ist, line, '\n')) {

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


                dataPositions.insert(std::make_pair(dataName, ist.tellg()));
                //ist.seekg(ist.tellg());
            }
        }
        ist.clear();
        return dataPositions;
    }

    template<typename T>
    static void readData(std::istream& ist, DataContainer<T, MeshDimension>& data) {

        std::map<std::string, std::istream::pos_type> dataPositions = indexData(ist);

        readCellData<typename DefaultIOTraits<T>::traitsType>::read(ist, data, dataPositions);
    }
};


#endif // VTKMESHDATAREADER_H
