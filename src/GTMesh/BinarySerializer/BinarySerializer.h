#ifndef BINARYSERIALIZER_H
#define BINARYSERIALIZER_H
#include <vector>
#include <cstdint>
#include <fstream>

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

    template<typename T, typename... TupleTraits>
    static void addNext(std::vector<Byte> &binaryDataContainer,
                        const T &data,
                        const std::tuple<TupleTraits...> &tupleTraits);

    template<typename T>
    static void readNext(std::vector<Byte>::const_iterator &binaryDataIterator, T &data);

    template<typename T, typename... TupleTraits>
    static void readNext(std::vector<Byte>::const_iterator &binaryDataIterator,
                         T &data,
                         const std::tuple<TupleTraits...> &tupleTraits);

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
        resetIteratorIfNotValid();
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

    /**
     * @brief This function clears the contained vector. Note that clear does not
     * free the memory. Thus, it is efficient to reuse serializer especially for sending
     * repetetive messages with the same data length.
     */
    void clear() {
        mData.clear();
        mDataIterator = {};
    }

    const ByteContainer& getData() {
        return mData;
    }

    void saveToFile(const std::string &fileName)
    {
        std::ofstream file(fileName, std::ofstream::binary);
        if (file) {
            file.write(reinterpret_cast<char *>(mData.data()), mData.size());
            file.close();
        } else {
            throw std::runtime_error("unable to open file " + fileName);
        }
    }

    void loadFromFile(const std::string &fileName) {
        std::ifstream file(fileName, std::ifstream::binary);
        if (file) {
            mData = ByteContainer(std::istreambuf_iterator<char>(file),
                                  std::istreambuf_iterator<char>());
        }
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

template<typename VarType, typename TupleTraits, typename Serializer, typename = void>
struct IsSerializableByWithBindedTraits : public std::false_type
{};

template<typename VarType, typename TupleTraits, typename Serializer>
struct IsSerializableByWithBindedTraits<
    VarType,
    TupleTraits,
    Serializer,
    decltype(Serializer::serialize(std::declval<BinarySerializer::ByteContainer &>(),
                                   std::declval<const VarType &>(),
                                   std::declval<const TupleTraits &>()))>
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

template<typename VarType, typename TupleTraits, typename Serializer, typename = void>
struct IsDeserializableByWithBindedTraits : public std::false_type
{};

template<typename VarType, typename TupleTraits, typename Serializer>
struct IsDeserializableByWithBindedTraits<
    VarType,
    TupleTraits,
    Serializer,
    decltype(Serializer::deserialize(std::declval<std::vector<unsigned char>::const_iterator &>(),
                                     std::declval<VarType &>(),
                                     std::declval<const TupleTraits &>()))> : public std::true_type
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

template<typename TupleType, typename Serializer>
struct IsSerializableByWithTraits
{};


template<typename VarType, typename TraitsTuple, typename Serializer>
struct IsSerializableByWithTraits<std::tuple<VarType, TraitsTuple>, Serializer>
    : Impl::IsSerializableByWithBindedTraits<VarType, TraitsTuple, Serializer>
{};


template<typename VarType, typename TraitsTuple, typename... Serializers>
using SelectSerializerWithTraits = ClassSelector<std::tuple<VarType, TraitsTuple>, IsSerializableByWithTraits, Serializers...>;

template<typename T, typename ... TupleTraits>
void BinarySerializer::addNext(std::vector<BinarySerializer::Byte> &binaryDataContainer, const T &data,
                               const std::tuple<TupleTraits...>& traits)
{
    SelectSerializerWithTraits<T, std::tuple<TupleTraits...>,
                               SerializeCustom,
                               SerializeSimple,
                               SerializeCompactContainer,
                               SerializeTraitedClass,
                               SerializeIterable,
                               SerializeIndexable,
                               SerializeTuple>::SelectedClass::serialize(binaryDataContainer,
                                                                         data,
                                                                         traits);
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

template<typename TupleType, typename Serializer>
struct IsDeserializableByWithTraits
{};


template<typename VarType, typename TraitsTuple, typename Serializer>
struct IsDeserializableByWithTraits<std::tuple<VarType, TraitsTuple>, Serializer>
    : Impl::IsDeserializableByWithBindedTraits<VarType, TraitsTuple, Serializer>
{};


template<typename VarType, typename TraitsTuple, typename... Serializers>
using SelectDeserializerWithTraits = ClassSelector<std::tuple<VarType, TraitsTuple>, IsDeserializableByWithTraits, Serializers...>;

template<typename T, typename ... TupleTraits>
void BinarySerializer::readNext(std::vector<Byte>::const_iterator& binaryDataIterator, T& data,
                                const std::tuple<TupleTraits...>& traits)
{
    SelectDeserializerWithTraits<T, std::tuple<TupleTraits...>,
                                 SerializeCustom,
                                 SerializeSimple,
                                 SerializeCompactContainer,
                                 SerializeTraitedClass,
                                 DeserializeAsociativeMap,
                                 DeserializeAsociativeSet,
                                 SerializeIterable,
                                 SerializeIndexable,
                                 SerializeTuple>::SelectedClass::deserialize(binaryDataIterator,
                                                                             data,
                                                                             traits);
}

#endif // BINARYSERIALIZER_H

