#ifndef SERIALIZECUSTOM_H
#define SERIALIZECUSTOM_H
#include "../BinarySerializer.h"
#include <GTMesh/Traits/CustomTypeTraits.h>

struct SerializeCustom {

    template <typename T>
    using IsBitviseSerializable = Impl::void_t<decltype(BinarySerialize(std::declval<BinarySerializer::ByteContainer&>(), std::declval<const T&>()))>;


    template <typename T>
    using IsBitviseDeserializable = Impl::void_t<decltype(BinaryDeserialize(std::declval<BinarySerializer::ByteContainerIterator&>(), std::declval<T&>()))>;


    template <typename T>
    static
        IsBitviseSerializable<T>
        serialize(BinarySerializer::ByteContainer& dataContainer, const T& data, ...) {
        BinarySerialize(dataContainer, data);
    }

    template <typename T>
    static
        IsBitviseDeserializable<T>
        deserialize(BinarySerializer::ByteContainerIterator& dataIterator, T& data, ...) {
        BinaryDeserialize(dataIterator, data);
    }
};
#endif // SERIALIZECUSTOM_H
