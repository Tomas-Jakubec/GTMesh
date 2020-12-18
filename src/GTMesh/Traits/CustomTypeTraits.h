#ifndef CUSTOMTRAITS_H
#define CUSTOMTRAITS_H
#include "Traits.h"
#include <iostream>
#include <type_traits>

namespace Impl {
#ifdef _MSC_VER
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
#else
template<typename...>
using void_t = void;

template<typename...>
using bool_t = bool;
#endif
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
struct is_tnl_indexable : public std::false_type
{};

template<typename T1>
struct is_tnl_indexable<
    T1,
    void_t<decltype(std::declval<const T1 &>()[0]), decltype(std::declval<const T1 &>().getSize())>>
    : public std::true_type
{};

template<typename T1, typename VOID = void>
struct is_traits : public std::false_type
{};

template<typename T1>
struct is_traits<T1, void_t<decltype(T1::isTraits)>> : public std::true_type
{};

template<typename T1, typename VOID = void>
struct has_default_traits : public std::false_type
{};

template<typename T1>
struct has_default_traits<T1, void_t<decltype(DefaultTraits<T1>::getTraits())>>
    : public std::true_type
{};

template<typename T1, typename = void>
struct has_default_io_traits : public std::false_type
{};

template<typename T1>
struct has_default_io_traits<T1, void_t<decltype(DefaultIOTraits<T1>::getTraits())>>
    : public std::true_type
{};

template<typename T1, typename VOID = void>
struct has_default_arithmetic_traits : public std::false_type
{};

template<typename T1>
struct has_default_arithmetic_traits<T1, void_t<decltype(DefaultArithmeticTraits<T1>::getTraits())>>
    : public std::true_type
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
struct IsTNLIndexable : public Impl::is_tnl_indexable<T1>
{};

template<typename T>
struct IsTraits : public Impl::is_traits<T>
{};

template<typename T>
struct HasDefaultTraits : public Impl::has_default_traits<T>
{};

template<typename T>
struct HasDefaultIOTraits : public Impl::has_default_io_traits<T>
{};

template<typename T>
struct HasDefaultArithmeticTraits : public Impl::has_default_arithmetic_traits<T>
{};

#endif // CUSTOMTRAITS_H
