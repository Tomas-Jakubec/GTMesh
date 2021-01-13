#ifndef SERIALIZEINDEXABLE_H
#define SERIALIZEINDEXABLE_H
#include <GTMesh/Traits/CustomTypeTraits.h>
#include "../BinarySerializer.h"
struct SerializeIndexable {
    template <typename T, typename ..., std::enable_if_t<IsIndexable<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data){
        using size_type = decltype(data.size());

        BinarySerializer::addNext(dataContainer, data.size());
        size_type size = data.size();

        for (size_type i = 0; i < size; i++) {
            BinarySerializer::addNext(dataContainer, data[i]);
        }
    }



    template <typename T, typename ..., std::enable_if_t<IsIndexable<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){
        using size_type = decltype(data.size());
        size_type size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        if (data.size() != size) {
            data.resize(size);
        }

        for (size_type i = 0; i < size; i++) {
            BinarySerializer::readNext(dataIterator, data[i]);
        }
    }

};
#endif // SERIALIZEINDEXABLE_H
