#ifndef CUSTOMTRAITS_H
#define CUSTOMTRAITS_H
#include <iostream>
#include <type_traits>
#include "Traits.h"

namespace Impl {



template <typename T1, typename T2 = void>
struct __is_exportable : public std::false_type {

};

template <typename T1>
struct __is_exportable<T1, typename std::enable_if<std::is_class<
        typename std::remove_reference<decltype(std::cerr << std::declval<const T1&>())>::type
        >::value>::type> : public std::true_type {};







template <typename T1, typename T2 = void>
struct __is_iterable : public std::false_type {};

template <typename T1>
struct __is_iterable<T1, typename std::enable_if<
            !std::is_same<
                decltype (std::declval<const T1&>().begin()),
                void
            >::value &&
            !std::is_same<
                decltype (std::declval<const T1&>().end()),
                void
            >::value
        >::type> : public std::true_type {};


template <typename T1, typename T2 = void>
struct __is_indexable : public std::false_type {};

template <typename T1>
struct __is_indexable<T1, typename std::enable_if<
            !std::is_same<
                decltype (std::declval<const T1&>()[0]),
                void
            >::value &&
            !std::is_same<
                decltype (std::declval<const T1&>().size()),
                void
            >::value
    >::type> : public std::true_type {};


template <typename T1, typename T2 = void>
struct __is_tnl_indexable : public std::false_type {};

template <typename T1>
struct __is_tnl_indexable<T1, typename std::enable_if<
            !std::is_same<
                decltype (std::declval<const T1&>()[0]), void
            >::value &&
            !std::is_same<
                decltype (std::declval<const T1&>().getSize()), void
            >::value
    >::type> : public std::true_type {};


template <typename T1, typename VOID = void>
struct __has_default_traits : public std::false_type {};


template <typename T1>
struct __has_default_traits<
        T1,
        typename std::enable_if<
            noexcept(Traits<T1>::getTraits)
        >::type
        > : public std::true_type {};



template <typename T1, typename VOID = void>
struct __has_default_io_traits : public std::false_type {};


template <typename T1>
struct __has_default_io_traits<
        T1,
        typename std::enable_if<
            noexcept(DefaultIOTraits<T1>::getTraits)
        >::type
        > : public std::true_type {};



template <typename T1, typename VOID = void>
struct __has_default_arithmetic_traits : public std::false_type {};


template <typename T1>
struct __has_default_arithmetic_traits<
        T1,
        typename std::enable_if<
            noexcept(DefaultArithmeticTraits<T1>::getTraits)
        >::type
        > : public std::true_type {};

} // Impl namespace


template <typename T1>
struct IsExportable : public Impl::__is_exportable<T1>
{};


template <typename T1>
struct IsIterable : public Impl::__is_iterable<T1>
{};

template <typename T1>
struct IsIndexable : public Impl::__is_indexable<T1>
{};

template <typename T1>
struct IsTNLIndexable : public Impl::__is_tnl_indexable<T1>
{};

template<typename T>
struct HasDefaultTraits : public Impl::__has_default_traits<T> {};

template<typename T>
struct HasDefaultIOTraits : public Impl::__has_default_io_traits<T> {};

template<typename T>
struct HasDefaultArithmeticTraits : public Impl::__has_default_arithmetic_traits<T> {};

#endif // CUSTOMTRAITS_H
