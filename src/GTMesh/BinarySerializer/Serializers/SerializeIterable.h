#ifndef SERIALIZEITERABLE_H
#define SERIALIZEITERABLE_H
#include <GTMesh/Traits/CustomTypeTraits.h>
#include "../BinarySerializer.h"

struct SerializeIterable {
    template <typename T, typename ..., std::enable_if_t<IsIterable<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data){
        BinarySerializer::addNext(dataContainer, data.size());
        for (const auto& val : data) {
            BinarySerializer::addNext(dataContainer, val);
        }
    }



    template <typename T, typename ..., std::enable_if_t<IsIterable<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){

        decltype(data.size()) size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        if (data.size() != size) {
            data.resize(size);
        }

        for (auto & val : data) {
            BinarySerializer::readNext(dataIterator, val);
        }
    }

    // TODO think how to deserialize lists and maps
};

struct DeserializeAsociativeMap {
    template <typename T, typename ..., std::enable_if_t<IsIterable<T>::value, bool> = true>
    static Impl::void_t<typename T::key_type,typename T::mapped_type>
        deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){

        decltype(data.size()) size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        for (decltype(data.size()) i = 0; i < size; i ++) {
            std::pair<typename T::key_type,typename T::mapped_type> val;
            BinarySerializer::readNext(dataIterator, val);
            data.insert(val);
        }
    }
};

struct DeserializeAsociativeSet {
    template <typename T, typename ..., std::enable_if_t<IsIterable<T>::value, bool> = true>
    static Impl::void_t<typename T::key_type>
    deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){

        decltype(data.size()) size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        for (decltype(data.size()) i = 0; i < size; i ++) {
            typename T::key_type val;
            BinarySerializer::readNext(dataIterator, val);
            data.insert(val);
        }
    }
};

#endif // SERIALIZEITERABLE_H
