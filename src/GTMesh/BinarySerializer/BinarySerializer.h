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

    template<typename T, typename ... Rest>
    void write(const T& data, const Rest&... additionalData) {
        write(data);
        write(additionalData...);
    }

    template<typename T>
    void write(const T& data) {
        addNext(mData, data);
    }

    template <typename T>
    BinarySerializer& operator<<(const T& value) {
        this->write(value);
        return *this;
    }

    template <typename T, typename ... Rest>
    void read(T& data, Rest&... additionalData) {
        read(data);
        read(additionalData...);
    }

    template <typename T>
    void read(T& data) {
        resetIteratorIfNotValid();
        readNext(mDataIterator, data);
    }

    template <typename T>
    T read() {
        return readNext<T>(mDataIterator);
    }

    template <typename T>
    BinarySerializer& operator>>(T& value) {
        this->read(value);
        return *this;
    }

    void resetIterator() {
        mDataIterator = mData.cbegin();
    }

    bool isIteratorValid() {
        return mData.begin() <= mDataIterator && mData.end() > mDataIterator;
    }

    void resetIteratorIfNotValid() {
        if (!isIteratorValid()) {
            resetIterator();
        }
    }

    void clear() {
        mData.clear();
        mDataIterator = {};
    }

    const ByteContainer& getData() {
        return mData;
    }
private:
    ByteContainer mData;
    ByteContainerIterator mDataIterator;
};

#include "Serializers/SerializeCustom.h"
#include "Serializers/SerializeSimple.h"
#include "Serializers/SerializeTraitedClass.h"
#include "Serializers/SerializeIterable.h"
#include "Serializers/SerializeIndexable.h"
#include "Serializers/SerializeTuple.h"
#include "Serializers/SerializeCompactContainer.h"
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
                                  std::declval<VarType &>()))> : public std::true_type
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
    SelectSerializer<T,
                     SerializeCustom,
                     SerializeSimple,
                     SerializeCompactContainer,
                     SerializeTraitedClass,
                     SerializeIterable,
                     SerializeIndexable,
                     SerializeTuple>::SelectedClass::serialize(binaryDataContainer, data);
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
                       SerializeCompactContainer,
                       SerializeTraitedClass,
                       DeserializeAsociativeMap,
                       DeserializeAsociativeSet,
                       SerializeIterable,
                       SerializeIndexable,
                       SerializeTuple>::SelectedClass::deserialize(binaryDataIterator, data);
}

#endif // BINARYSERIALIZER_H

