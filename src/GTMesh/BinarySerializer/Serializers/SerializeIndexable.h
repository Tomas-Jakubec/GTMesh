#ifndef SERIALIZEINDEXABLE_H
#define SERIALIZEINDEXABLE_H
#include <GTMesh/Traits/CustomTypeTraits.h>
#include "../BinarySerializer.h"

namespace Interface {
/**
 * Generic template for indexable interface
 */
template <typename T, typename = void>
struct Indexable : public std::false_type {};

/**
 * Indexable interface for stl containers.
 */
template<typename T>
struct Indexable<T, std::enable_if_t<::IsIndexable<T>::value>> : public std::true_type
{
    using size_type = decltype(std::declval<const T &>().size());

    static size_type size(const T &val) { return val.size(); }

    template<typename _T = T, std::enable_if_t<IsResizable<_T>::value, bool> = true>
    static void resize(T &val, size_type new_size)
    {
        val.resize(new_size);
    }

    /**
     * Some types like std::array do not have the resize member function.
     * However, it still has subscript operator.
     */
    template<typename _T = T, std::enable_if_t<!IsResizable<_T>::value, bool> = true>
    static void resize(T &/*val*/, size_type /*new_size*/)
    {}
};

/**
 * Indexable interface for TNL containers.
 */
template<typename T>
struct Indexable<T, std::enable_if_t<::IsTNLIndexable<T>::value>> : public std::true_type
{
    using size_type = decltype(std::declval<const T &>().getSize());

    static size_type size(const T &val) { return val.getSize(); }

    template<typename _T = T, std::enable_if_t<IsTNLResizable<_T>::value, bool> = true>
    static void resize(T &val, size_type new_size)
    {
        val.setSize(new_size);
    }

    template<typename _T = T, std::enable_if_t<!IsTNLResizable<_T>::value, bool> = true>
    static void resize(T &/*val*/, size_type /*new_size*/)
    {}
};
}

struct SerializeIndexable {
    template <typename T, typename ..., std::enable_if_t<Interface::Indexable<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data){

        auto size = Interface::Indexable<T>::size(data);
        BinarySerializer::addNext(dataContainer, size);

        for (decltype(size) i = 0; i < size; i++) {
            BinarySerializer::addNext(dataContainer, data[i]);
        }
    }

    template <typename T, typename ... TupleTraits, std::enable_if_t<Interface::Indexable<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data,
                          const std::tuple<TupleTraits...>& traits){

        auto size = Interface::Indexable<T>::size(data);
        BinarySerializer::addNext(dataContainer, size);

        for (decltype(size) i = 0; i < size; i++) {
            BinarySerializer::addNext(dataContainer, data[i], traits);
        }
    }

    template <typename T, typename ..., std::enable_if_t<IsIndexable<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){
        using size_type = decltype(data.size());
        size_type size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        if (Interface::Indexable<T>::size(data) != size) {
            Interface::Indexable<T>::resize(data, size);
        }

        for (size_type i = 0; i < size; i++) {
            BinarySerializer::readNext(dataIterator, data[i]);
        }
    }


    template <typename T, typename ... TupleTraits, std::enable_if_t<IsIndexable<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data,
                            const std::tuple<TupleTraits...>& traits){
        using size_type = decltype(data.size());
        size_type size = BinarySerializer::readNext<decltype(data.size())>(dataIterator);

        if (Interface::Indexable<T>::size(data) != size) {
            Interface::Indexable<T>::resize(data, size);
        }

        for (size_type i = 0; i < size; i++) {
            BinarySerializer::readNext(dataIterator, data[i], traits);
        }
    }

};

#endif // SERIALIZEINDEXABLE_H
