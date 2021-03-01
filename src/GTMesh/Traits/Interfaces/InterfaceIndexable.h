#ifndef INTERFACEINDEXABLE_H
#define INTERFACEINDEXABLE_H
#include <type_traits>
#include "../CustomTypeTraits.h"

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
    using sizeType = decltype(std::declval<const T &>().size());

    using valueType = typename std::decay<decltype(std::declval<const T &>()[0])>::type;

    static sizeType size(const T &val) { return val.size(); }

    static const valueType& getElement(const T& array, const sizeType& index) {
        return array[index];
    }

    static void setElement(T& array, const sizeType& index, const valueType& newValue) {
        return array[index] = newValue;
    }

    static valueType& getElementRef(T& array, const sizeType& index) {
        return array[index];
    }

    template<typename _T = T, std::enable_if_t<IsResizable<_T>::value, bool> = true>
    static void resize(T &val, sizeType new_size)
    {
        val.resize(new_size);
    }

    /**
     * Some types like std::array do not have the resize member function.
     * However, it still has subscript operator.
     */
    template<typename _T = T, std::enable_if_t<!IsResizable<_T>::value, bool> = true>
    static void resize(T &/*val*/, sizeType /*new_size*/)
    {}
};

/**
 * Indexable interface for TNL containers.
 */
template<typename T>
struct Indexable<T, std::enable_if_t<::IsTNLIndexable<T>::value>> : public std::true_type
{
    using sizeType = decltype(std::declval<const T &>().getSize());
    using valueType = typename std::decay<decltype(std::declval<const T &>()[0])>::type;

    static sizeType size(const T &val) { return val.getSize(); }

    static const valueType& getElement(const T& array, const sizeType& index) {
        return array[index];
    }

    static void setElement(T& array, const sizeType& index, const valueType& newValue) {
        return array[index] = newValue;
    }

    static valueType& getElementRef(T& array, const sizeType& index) {
        return array[index];
    }

    template<typename _T = T, std::enable_if_t<IsTNLResizable<_T>::value, bool> = true>
    static void resize(T &val, sizeType new_size)
    {
        val.setSize(new_size);
    }

    template<typename _T = T, std::enable_if_t<!IsTNLResizable<_T>::value, bool> = true>
    static void resize(T &/*val*/, sizeType /*new_size*/)
    {}
};
}
#endif // INTERFACEINDEXABLE_H
