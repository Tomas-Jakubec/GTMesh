#ifndef CUSTOMTRAITS_H
#define CUSTOMTRAITS_H
#include <iostream>
#include <type_traits>

namespace Impl {
template<typename... T> struct make_void {
    using type = void;
};

template<typename... T>
using void_t = typename make_void<T...>::type;

template<typename... T> struct make_bool {
    using type = bool;
};

template<typename... T>
using bool_t = typename make_bool<T...>::type;

template<typename T1, typename T2 = void>
struct is_exportable : public std::false_type
{};

template<typename T1>
struct is_exportable<T1, void_t<decltype(std::cerr << std::declval<const T1 &>())>>
    : public std::true_type
{};

template<typename T1, typename T2 = void>
struct is_iterable : public std::false_type
{};

template<typename T1>
struct is_iterable<
    T1,
    void_t<decltype(std::declval<const T1 &>().begin()), decltype(std::declval<const T1 &>().end())>>
    : public std::true_type
{};

template<typename T1, typename T2 = void>
struct is_indexable : public std::false_type
{};

template<typename T1>
struct is_indexable<
    T1,
    void_t<decltype(std::declval<const T1 &>()[0]), decltype(std::declval<const T1 &>().size())>>
    : public std::true_type
{};

template<typename T1, typename T2 = void>
struct is_resizable : public std::false_type
{};

template<typename T1>
struct is_resizable<
    T1,
    void_t<decltype(std::declval<T1 &>().resize(0))>>
    : public std::true_type
{};

template<typename T1, typename T2 = void>
struct is_tnl_indexable : public std::false_type
{};

template<typename T1>
struct is_tnl_indexable<
    T1,
    void_t<decltype(std::declval<const T1 &>()[0]), decltype(std::declval<const T1 &>().getSize())>>
    : public std::true_type
{};

template<typename T1, typename T2 = void>
struct is_tnl_resizable : public std::false_type
{};

template<typename T1>
struct is_tnl_resizable<
    T1,
    void_t<decltype(std::declval<T1 &>().setSize(0))>>
    : public std::true_type
{};

template<typename T1, typename VOID = void>
struct is_traits : public std::false_type
{};

template<typename T1>
struct is_traits<T1, void_t<decltype(T1::isTraits)>> : public std::true_type
{};



} // namespace Impl

template<typename T1>
struct IsExportable : public Impl::is_exportable<T1>
{};

template<typename T1>
struct IsIterable : public Impl::is_iterable<T1>
{};

template<typename T1>
struct IsIndexable : public Impl::is_indexable<T1>
{};

template<typename T1>
struct IsResizable : public Impl::is_resizable<T1>
{};

template<typename T1>
struct IsTNLIndexable : public Impl::is_tnl_indexable<T1>
{};

template<typename T1>
struct IsTNLResizable : public Impl::is_tnl_resizable<T1>
{};

template<typename T>
struct IsTraits : public Impl::is_traits<T>
{};

#endif // CUSTOMTRAITS_H
