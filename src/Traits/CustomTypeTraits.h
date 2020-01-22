#ifndef CUSTOMTRAITS_H
#define CUSTOMTRAITS_H
#include <iostream>
#include <type_traits>
#include "Traits.h"

namespace Detail {



template <typename T1, typename T2 = void>
struct __is_exportable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_exportable<T1, typename std::enable_if<std::is_class<
        typename std::remove_reference<decltype(std::cerr << std::declval<const T1&>())>::type
        >::value>::type> : public std::integral_constant<bool, true> {

};





template <typename T1, typename T2 = void>
struct __is_iterable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_iterable<T1, typename std::enable_if<!std::is_same<
        decltype(std::declval<const T1&>().begin()),
        void
     >::value>::type> : public std::integral_constant<bool, true> {

};


template <typename T1, typename T2 = void>
struct __is_indexable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_indexable<T1, typename std::enable_if<
            !std::is_same<
                decltype (std::declval<const T1&>()[0]), void
            >::value &&
            !std::is_same<
                decltype (std::declval<const T1&>().size()), void
            >::value
    >::type> : public std::integral_constant<bool, true> {

};


template <typename T1, typename T2 = void>
struct __is_tnl_indexable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_tnl_indexable<T1, typename std::enable_if<
            !std::is_same<
                decltype (std::declval<const T1&>()[0]), void
            >::value &&
            !std::is_same<
                decltype (std::declval<const T1&>().getSize()), void
            >::value
    >::type> : public std::integral_constant<bool, true> {

};


template <typename T1, typename VOID = void>
struct __has_default_traits : public std::integral_constant<bool, false> {};


template <typename T1>
struct __has_default_traits<
        T1,
        typename std::enable_if<
            Traits<T1>::is_specialized
        >::type
        > : public std::integral_constant<bool, true> {};



template <typename T1, typename VOID = void>
struct __has_default_io_traits : public std::integral_constant<bool, false> {};


template <typename T1>
struct __has_default_io_traits<
        T1,
        typename std::enable_if<
            DefaultIOTraits<T1>::is_specialized
        >::type
        > : public std::integral_constant<bool, true> {};



template <typename T1, typename VOID = void>
struct __has_default_arithmetic_traits : public std::integral_constant<bool, false> {};


template <typename T1>
struct __has_default_arithmetic_traits<
        T1,
        typename std::enable_if<
            DefaultArithmeticTraits<T1>::is_specialized
        >::type
        > : public std::integral_constant<bool, true> {};

}


template <typename T1>
struct IsExportable : public Detail::__is_exportable<T1>
{};


template <typename T1>
struct IsIterable : public Detail::__is_iterable<T1>
{};

template <typename T1>
struct IsIndexable : public Detail::__is_indexable<T1>
{};

template <typename T1>
struct IsTNLIndexable : public Detail::__is_tnl_indexable<T1>
{};

template<typename T>
struct HasDefaultTraits : public Detail::__has_default_traits<T> {};

template<typename T>
struct HasDefaultIOTraits : public Detail::__has_default_io_traits<T> {};

template<typename T>
struct HasDefaultArithmeticTraits : public Detail::__has_default_arithmetic_traits<T> {};

#endif // CUSTOMTRAITS_H
