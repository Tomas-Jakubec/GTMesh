#ifndef BINARYSERIALIZER_H
#define BINARYSERIALIZER_H
#include <vector>
#include <cstdint>

// TODO možná změnit jméno
/**
 * @brief The BinarySerializer class holds binary array and is able to
 * store binary data as well as serialize and deserialize data
 */
class BinarySerializer {
public:
    using Byte = uint8_t;
    using ByteContainer = std::vector<Byte>;
    using ByteContainerIterator = ByteContainer::const_iterator;


    template<typename T>
    static void addNext(std::vector<Byte>& binaryDataContainer, const T& data);

    template <typename T>
    static void readNext(std::vector<Byte>::const_iterator& binaryDataIterator, T& data);

    template <typename T>
    static T readNext(std::vector<Byte>::const_iterator& binaryDataIterator) {
        T data;
        readNext(binaryDataIterator, data);
        return data;
    }

    template<typename T>
    void maddNext(const T& data) {
        addNext(mData, data);
    }

    template <typename T>
    void mreadNext(T& data) {
        readNext(mDataIterator, data);
    }

    template <typename T>
    T mreadNext() {
        return readNext<T>(mDataIterator);
    }

    void resetIterator() {
        mDataIterator = mData.cbegin();
    }
//private:
    ByteContainer mData;
    ByteContainerIterator mDataIterator;
};

#include "Serializers/SerializeCustom.h"
#include "Serializers/SerializeSimple.h"
#include "Serializers/SerializeIterable.h"
#include "Serializers/SerializeTuple.h"
#include <GTMesh/Utils/ClassSelector.h>

namespace Impl {
template<typename VarType, typename Serializer, typename = void>
struct IsSerializableBy : public std::false_type
{};

template<typename VarType, typename Serializer>
struct IsSerializableBy<VarType,
                        Serializer,
                        decltype(Serializer::serialize(std::declval<std::vector<unsigned char> &>(),
                                                    std::declval<const VarType &>()))>
    : public std::true_type
{};

template<typename VarType, typename Serializer, typename = void>
struct IsDeserializableBy : public std::false_type
{};

template<typename VarType, typename Serializer>
struct IsDeserializableBy<
    VarType,
    Serializer,
    decltype(Serializer::deserialize(std::declval<std::vector<unsigned char>::const_iterator &>(),
                                  std::declval<const VarType &>()))> : public std::true_type
{};

} // Impl namespace

template<typename VarType, typename Serializer>
struct IsSerializableBy: public Impl::IsSerializableBy<VarType, Serializer>
{};

template<typename VarType, typename... Serializers>
using SelectSerializer = ClassSelector<VarType, IsSerializableBy, Serializers...>;

template<typename T>
void BinarySerializer::addNext(std::vector<BinarySerializer::Byte> &binaryDataContainer, const T &data)
{
    SelectSerializer<T, SerializeCustom, SerializeSimple, SerializeIterable, SerializeTuple>::
        SelectedClass::serialize(binaryDataContainer, data);
}

template<typename VarType, typename Serializer>
struct IsDeserializableBy: public Impl::IsDeserializableBy<VarType, Serializer>
{};

template<typename VarType, typename... Serializers>
using SelectDeserializer = ClassSelector<VarType, IsDeserializableBy, Serializers...>;

template<typename T>
void BinarySerializer::readNext(std::vector<Byte>::const_iterator& binaryDataIterator, T& data)
{
    SelectDeserializer<T,
                       SerializeCustom,
                       SerializeSimple,
                       DeserializeAsociativeMap,
                       DeserializeAsociativeSet,
                       SerializeIterable,
                       SerializeTuple>::SelectedClass::deserialize(binaryDataIterator, data);
}

#endif // BINARYSERIALIZER_H

