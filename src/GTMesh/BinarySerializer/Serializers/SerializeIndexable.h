#ifndef SERIALIZEINDEXABLE_H
#define SERIALIZEINDEXABLE_H
#include <GTMesh/Traits/CustomTypeTraits.h>
#include <GTMesh/Traits/Interfaces/InterfaceIndexable.h>
#include "../BinarySerializer.h"


struct SerializeIndexable {
    template <typename T, typename ..., std::enable_if_t<Interface::Indexable<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data){

        auto size = Interface::Indexable<T>::size(data);
        BinarySerializer::addNext(dataContainer, size);

        for (decltype(size) i = 0; i < size; i++) {
            BinarySerializer::addNext(dataContainer, Interface::Indexable<T>::getElement(data, i));
        }
    }

    template <typename T, typename ... TupleTraits, std::enable_if_t<Interface::Indexable<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data,
                          const std::tuple<TupleTraits...>& traits){

        auto size = Interface::Indexable<T>::size(data);
        BinarySerializer::addNext(dataContainer, size);

        for (decltype(size) i = 0; i < size; i++) {
            BinarySerializer::addNext(dataContainer, Interface::Indexable<T>::getElement(data, i), traits);
        }
    }

    template <typename T, typename ..., std::enable_if_t<Interface::Indexable<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){
        using size_type = decltype(data.size());
        size_type size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        if (Interface::Indexable<T>::size(data) != size) {
            Interface::Indexable<T>::resize(data, size);
        }

        for (size_type i = 0; i < size; i++) {
            BinarySerializer::readNext(dataIterator, Interface::Indexable<T>::getElementRef(data, i));
        }
    }


    template <typename T, typename ... TupleTraits, std::enable_if_t<Interface::Indexable<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data,
                            const std::tuple<TupleTraits...>& traits){
        using size_type = decltype(data.size());
        size_type size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        if (Interface::Indexable<T>::size(data) != size) {
            Interface::Indexable<T>::resize(data, size);
        }

        for (size_type i = 0; i < size; i++) {
            BinarySerializer::readNext(dataIterator, Interface::Indexable<T>::getElementRef(data, i), traits);
        }
    }

};

#endif // SERIALIZEINDEXABLE_H
