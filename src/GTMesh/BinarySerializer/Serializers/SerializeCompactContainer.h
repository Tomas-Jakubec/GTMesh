#ifndef SERIALIZECOMPACTCONTAINER_H
#define SERIALIZECOMPACTCONTAINER_H
#include <GTMesh/Traits/CustomTypeTraits.h>
#include "../BinarySerializer.h"
#include "SerializeSimple.h"
#include <GTMesh/Debug/Debug.h>
namespace Interface {
/**
 * Generic template for indexable interface
 */
template<typename T, typename = void>
struct Compact : public std::false_type
{};

template<typename T>
struct Compact<T, std::enable_if_t<Impl::IsSimpleSerializable_v<T>>> : public std::true_type
{};

/**
 * Indexable interface for stl containers.
 */
template<typename T>
struct Compact<
    T,
    Impl::void_t<
        std::enable_if_t<::IsIndexable<T>::value
                         && Compact<std::decay_t<decltype(std::declval<const T &>()[0])>>::value>,
        decltype(std::declval<T &>().data())>> : public std::true_type
{
    using size_type = decltype(std::declval<const T &>().size());
    using value_type = std::decay_t<decltype(std::declval<const T &>()[0])>;

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
    static void resize(T & /*val*/, size_type /*new_size*/)
    {}

    static const value_type *data(const T &val) { return val.data(); }

    static value_type *data(T &val) { return val.data(); }
};

/**
 * Indexable interface for TNL containers.
 */
template<typename T>
struct Compact<
    T,
    void_t<std::enable_if_t<::IsTNLIndexable<T>::value
                            && Compact<std::decay_t<decltype(std::declval<const T &>()[0])>>::value>,
           decltype(std::declval<T &>().getData())>> : public std::true_type
{
    using size_type = decltype(std::declval<const T &>().getSize());
    using value_type = std::decay_t<decltype(std::declval<const T &>()[0])>;

    static size_type size(const T &val) { return val.getSize(); }

    template<typename _T = T, std::enable_if_t<IsTNLResizable<_T>::value, bool> = true>
    static void resize(T &val, size_type new_size)
    {
        val.setSize(new_size);
    }

    template<typename _T = T, std::enable_if_t<!IsTNLResizable<_T>::value, bool> = true>
    static void resize(T & /*val*/, size_type /*new_size*/)
    {}

    static const value_type *data(const T &val) { return val.getData(); }

    static value_type *data(T &val) { return val.getData(); }
};
} // namespace Interface

/**
 * @brief Accelerated serialization of containers which are aligned in memory similarly to std::vector.
 * The condition is that the contained type must be also simple binary serializable.
 */
struct SerializeCompactContainer {
    template <typename T, typename ..., std::enable_if_t<Interface::Indexable<T>::value && Interface::Compact<T>::value, bool> = true>
    static void serialize(std::vector<unsigned char>& dataContainer, const T& data){

        auto size = Interface::Compact<T>::size(data);
        /*
        dataContainer.resize(dataContainer.size() + (size * sizeof (Interface::Compact<T>::value_type)));
        */
        BinarySerializer::addNext(dataContainer, size);
        dataContainer.insert(dataContainer.end(),
                             reinterpret_cast<const BinarySerializer::Byte*>(Interface::Compact<T>::data(data)),
                             reinterpret_cast<const BinarySerializer::Byte*>(Interface::Compact<T>::data(data) + size));
    }

    template <typename T, typename ..., std::enable_if_t<Interface::Indexable<T>::value && Interface::Compact<T>::value, bool> = true>
    static void deserialize(std::vector<unsigned char>::const_iterator& dataIterator, T& data){
        auto size = Interface::Compact<T>::size(data);

        BinarySerializer::readNext(dataIterator, size);

        if (Interface::Compact<T>::size(data) != size) {
            Interface::Compact<T>::resize(data, size);
        }

        DBGVAR(size, Interface::Compact<T>::size(data));
        memcpy(Interface::Compact<T>::data(data), dataIterator.base(), size * sizeof (typename Interface::Compact<T>::value_type));
        dataIterator += size * sizeof (typename Interface::Compact<T>::value_type);
    }

};
#endif // SERIALIZECOMPACTCONTAINER_H
