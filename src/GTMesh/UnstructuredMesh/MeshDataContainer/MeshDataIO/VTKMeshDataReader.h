#ifndef VTKMESHDATAREADER_H
#define VTKMESHDATAREADER_H
#include "../../../Traits/Traits.h"
#include "../../../Traits/CustomTypeTraits.h"
#include "../MeshDataContainer.h"
#include "../../../Debug/Debug.h"
#include <istream>
#include <map>
#include <sstream>
#include <GTMesh/Debug/Printers/PrintTraitedClass.h>

template <unsigned int MeshDimension, typename IndexType>
class VTKMeshDataReader {

    static_assert (MeshDimension == 2 || MeshDimension == 3, "The VTK file format can represent data only in 2D or 3D");

    /**
     * @brief readColumn
     * reads a single column of traited data
     */
    template<typename T, unsigned int Index, typename ... TraitsArgs>
    static
        std::enable_if_t< IsIndexable<typename Traits<T, TraitsArgs...>::template type<Index>>::value &&
                         MeshDimension == 3 >
        readColumn( std::istream& ist,
                   DataContainer<T, MeshDimension> &data,
                   std::map<std::string, std::istream::pos_type>& dataPositions,
                   const Traits<T, TraitsArgs...>& traits)
    {

        ist.seekg(dataPositions[traits.template getName<Index>()]);

        typename Traits<T, TraitsArgs...>::template type<Index> value;

        for (IndexType i = 0; i < data.size(); i++) {
            for (unsigned int j = 0; j < traits.template getValue<Index>(data.at(i)).size(); j++){
                ist >> value[j];
            }
            traits.template setValue<Index>(data.at(i), value);
        }

    }


    template<typename T, unsigned int Index, typename ... TraitsArgs>
    static auto readColumn( std::istream& ist,
                           DataContainer<T, MeshDimension> &data,
                           std::map<std::string, std::istream::pos_type>& dataPositions,
                           const Traits<T, TraitsArgs...>& traits)
        -> typename std::enable_if<
            IsIndexable<typename DefaultIOTraits<T>::traitsType::template type<Index>>::value &&
            MeshDimension == 2
            >::type
    {

        ist.seekg(dataPositions[traits.template getName<Index>()]);

        typename Traits<T, TraitsArgs...>::template type<Index> value, dummy;

        for (IndexType i = 0; i < data.size(); i++) {
            for (unsigned int j = 0; j < traits.template getValue<Index>(data.at(i)).size(); j++){
                ist >> value[j];
            }

            ist >> dummy[0];

            traits.template setValue<Index>(data.at(i), value);
        }

    }


    template<typename T, unsigned int Index, typename ... TraitsArgs>
    static auto readColumn( std::istream& ist,
                           DataContainer<T, MeshDimension> &data,
                           std::map<std::string, std::istream::pos_type>& dataPositions,
                           const Traits<T, TraitsArgs...>& traits)
        -> typename std::enable_if<
            !IsIndexable<typename Traits<T, TraitsArgs...>::template type<Index>>::value
            >::type
    {


        ist.seekg(dataPositions[traits.template getName<Index>()]);

        typename Traits<T, TraitsArgs...>::template type<Index> value;

        for (IndexType i = 0; i < data.size(); i++){
            ist >> value;
            traits.template setValue<Index>(data.at(i), value);
        }

    }
private:
    // A helper class to be executed by the readFromStream member function of
    // MeshDataIterator
    struct ReadData {
        template<unsigned int Index,
                 typename T,
                 typename... TraitsArgs>
        static void exec(std::istream &ist,
                         DataContainer<T, MeshDimension> &data,
                         std::map<std::string, std::istream::pos_type> &dataPositions,
                         const Traits<T, TraitsArgs...> &traits = DefaultIOTraits<T>::getTraits())
        {
            readColumn<T, Index>(ist, data, dataPositions, traits);
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
    static void readFromStream(std::istream& ist, DataContainer<T, MeshDimension>& data) {

        std::map<std::string, std::istream::pos_type> dataPositions = indexData(ist);

        readData(ist, data, dataPositions);
    }
    // Search for importable containers

    template<unsigned int Index, bool OK = false>
    struct MeshDataIterator
    {
        template<typename T, unsigned int... Dimensions, typename... TraitsTuple>
        static typename std::enable_if<(!SelectTraits<typename MeshDataContainer<T, Dimensions...>::
                                                          template DataContainerType<Index>::type,
                                                      0,
                                                      TraitsTuple...>::valid)>::type
        readFromStream(std::istream &ist,
                       MeshDataContainer<T, Dimensions...> &data,
                       std::map<std::string, std::istream::pos_type> &dataPositions,
                       const std::tuple<TraitsTuple...> &tupTraits)
        {
            MeshDataIterator<Index - 1, OK>::readFromStream(ist, data, dataPositions, tupTraits);
        }

        template<typename T, unsigned int... Dimensions, typename... TraitsTuple>
        static typename std::enable_if<(SelectTraits<typename MeshDataContainer<T, Dimensions...>::
                                                         template DataContainerType<Index>::type,
                                                     0,
                                                     TraitsTuple...>::valid)>::type
        readFromStream(std::istream &ist,
                       MeshDataContainer<T, Dimensions...> &data,
                       std::map<std::string, std::istream::pos_type> &dataPositions,
                       const std::tuple<TraitsTuple...> &tupTraits)
        {
            using type =
                typename MeshDataContainer<T, Dimensions...>::template DataContainerType<Index>::type;

            MeshDataIterator<Index - 1, true>::writeToStream(ist, data, dataPositions, tupTraits);

            constexprFor<ReadData, SelectTraits<type, 0, TraitsTuple...>::TypeTraits::size()>(
                ist,
                data.template getDataByPos<Index>(),
                dataPositions,
                SelectTraits<type, 0, TraitsTuple...>::getTraitsInstance(tupTraits));
        }
    };

    template<bool OK>
    struct MeshDataIterator<0, OK>
    {
        template<typename T, unsigned int... Dimensions, typename... TraitsTuple>
        static typename std::enable_if<
            (!SelectTraits<
                typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type,
                0,
                TraitsTuple...>::valid)>::type
        readFromStream(std::istream &,
                       MeshDataContainer<T, Dimensions...> &,
                       std::map<std::string, std::istream::pos_type> &,
                       const std::tuple<TraitsTuple...> &)
        {
            static_assert(
                OK,
                "The mesh data container must have at least one DataContainer mapped to cells with "
                "traits for example using macro MAKE_ATTRIBUTE_TRAIT see header Traits.h");
        }

        template<typename T, unsigned int... Dimensions, typename... TraitsTuple>
        static typename std::enable_if<
            (SelectTraits<
                typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type,
                0,
                TraitsTuple...>::valid)>::type
        readFromStream(std::istream &ist,
                       MeshDataContainer<T, Dimensions...> &data,
                       std::map<std::string, std::istream::pos_type> &dataPositions,
                       const std::tuple<TraitsTuple...> &tupTraits)
        {
            using type =
                typename MeshDataContainer<T, Dimensions...>::template DataContainerType<0>::type;

            constexprFor<ReadData, SelectTraits<type, 0, TraitsTuple...>::TypeTraits::size()>(
                ist,
                data.template getDataByPos<0>(),
                dataPositions,
                SelectTraits<type, 0, TraitsTuple...>::getTraitsInstance(tupTraits));
        }
    };

public:
    template<typename T, unsigned int ...Dimensions, typename ... TraitsTypes>
    static void readFromStream( std::istream& ist,
                               MeshDataContainer<T, Dimensions...>& data,
                               const std::tuple<TraitsTypes...>& tupTraits = std::tuple<>()) {

        std::map<std::string, std::istream::pos_type> dataPositions = indexData(ist);

        MeshDataIterator<sizeof... (Dimensions) - 1>::readFromStream(ist, data, dataPositions, tupTraits);
    }
    template<typename T, unsigned int ...Dimensions, typename ... TraitsTypes>
    static void readFromStream( std::istream& ist,
                               const TraitsBinder<MeshDataContainer<T, Dimensions...>, TraitsTypes...>& data) {

        std::map<std::string, std::istream::pos_type> dataPositions = indexData(ist);

        MeshDataIterator<sizeof... (Dimensions) - 1>::readFromStream(ist, data.object, dataPositions, data.tupTraits);
    }
};


#endif // VTKMESHDATAREADER_H
