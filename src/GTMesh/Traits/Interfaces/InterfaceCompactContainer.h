#ifndef INTERFACECOMPACTCONTAINER_H
#define INTERFACECOMPACTCONTAINER_H
#include "InterfaceCompact.h"
#include "../CustomTypeTraits.h"
namespace Interface {

template<typename T, typename = void>
struct CompactContainer : public std::false_type
{};
/**
 * Indexable interface for stl containers.
 */
template<typename T>
struct CompactContainer<
    T,
    ::Impl::void_t<
        std::enable_if_t<::IsIndexable<T>::value
                         && Compact<std::decay_t<decltype(std::declval<const T &>()[0])>>::value>,
        decltype(std::declval<const T &>().data())>> : public std::true_type
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

    static const value_type *data(const T &val) { return &val[0]; }

    static value_type *data(T &val) { return &val[0]; }
};

/**
 * Indexable interface for TNL containers.
 */
template<typename T>
struct CompactContainer<
    T,
    ::Impl::void_t<std::enable_if_t<::IsTNLIndexable<T>::value
                            && Compact<std::decay_t<decltype(std::declval<const T &>()[0])>>::value>,
           decltype(std::declval<const T &>().getData())>> : public std::true_type
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

    static const value_type *data(const T &val) { return &val[0]; }

    static value_type *data(T &val) { return &val[0]; }
};
}

#endif // INTERFACECOMPACTCONTAINER_H
