#ifndef SERIALIZETUPLE_H
#define SERIALIZETUPLE_H
#include "../BinarySerializer.h"
#include <GTMesh/Utils/ConstexprFor.h>
#include <tuple>


struct SerializeTuple {
    /**
     * @brief Python-like print of dictionary pair of key and value.
     */
    template<typename T1, typename T2>
    static void serialize(BinarySerializer::ByteContainer& dataContainer, const std::pair<T1,T2>& data, ...)
    {
        BinarySerializer::addNext(dataContainer, data.first);
        BinarySerializer::addNext(dataContainer, data.second);
    }

    template<typename T1, typename T2>
    static void deserialize(BinarySerializer::ByteContainerIterator& dataIterator, std::pair<T1,T2>& data, ...)
    {
        BinarySerializer::readNext(dataIterator, data.first);
        BinarySerializer::readNext(dataIterator, data.second);
    }


    struct SerializeTupleHelper{
        template <unsigned int Index, typename ... TupleTypes>
        static void exec(BinarySerializer::ByteContainer& dataContainer, const std::tuple<TupleTypes...>& data) {
            BinarySerializer::addNext(dataContainer, std::get<Index>(data));
        }
    };

    template<typename ... TupleTypes>
    static void serialize(BinarySerializer::ByteContainer& dataContainer, const std::tuple<TupleTypes...>& data, ...)
    {
        constexprFor<SerializeTupleHelper, sizeof... (TupleTypes)>(dataContainer, data);
    }

    struct DeserializeTupleHelper{
        template <unsigned int Index, typename ... TupleTypes>
        static void exec(BinarySerializer::ByteContainerIterator& dataIterator, std::tuple<TupleTypes...>& data) {
            BinarySerializer::readNext(dataIterator, std::get<Index>(data));
        }
    };

    template<typename ... TupleTypes>
    static void deserialize(BinarySerializer::ByteContainerIterator& dataIterator, std::tuple<TupleTypes...>& data, ...)
    {
        constexprFor<DeserializeTupleHelper, sizeof... (TupleTypes)>(dataIterator, data);
    }
};
#endif // SERIALIZETUPLE_H
