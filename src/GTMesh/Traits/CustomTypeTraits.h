#ifndef CUSTOMTRAITS_H
#define CUSTOMTRAITS_H
#include <iostream>
#include <type_traits>
#include "Traits.h"

namespace Impl {

template<typename ...>
using void_t = void;

template<typename ...>
using bool_t = bool;


template <typename T1, typename T2 = void>
struct __is_exportable : public std::false_type {

};

template <typename T1>
struct __is_exportable< T1,
                        void_t<decltype(std::cerr << std::declval<const T1&>())> > : public std::true_type {};



template <typename T1, typename T2 = void>
struct __is_iterable : public std::false_type {};

template <typename T1>
struct __is_iterable< T1,
                      void_t<decltype(std::declval<const T1&>().begin()), decltype (std::declval<const T1&>().end())> > : public std::true_type {};


template <typename T1, typename T2 = void>
struct __is_indexable : public std::false_type {};

template <typename T1>
struct __is_indexable< T1,
                       void_t<decltype (std::declval<const T1&>()[0]), decltype (std::declval<const T1&>().size())> > : public std::true_type {};


template <typename T1, typename T2 = void>
struct __is_tnl_indexable : public std::false_type {};

template <typename T1>
struct __is_tnl_indexable< T1,
                           void_t<decltype (std::declval<const T1&>()[0]), decltype (std::declval<const T1&>().getSize())> > : public std::true_type {};



template <typename T1, typename VOID = void>
struct __is_traits : public std::false_type {};


template <typename T1>
struct __is_traits<
        T1,
        void_t<decltype(T1::isTraits)>
        > : public std::true_type {};


template <typename T1, typename VOID = void>
struct __has_default_traits : public std::false_type {};


template <typename T1>
struct __has_default_traits<
        T1,
        void_t<decltype(DefaultTraits<T1>::getTraits())>
        > : public std::true_type {};



template <typename T1, typename = void>
struct __has_default_io_traits : public std::false_type {};


template <typename T1>
struct __has_default_io_traits<
        T1,
        void_t<decltype(DefaultIOTraits<T1>::getTraits())>
        > : public std::true_type {};



template <typename T1, typename VOID = void>
struct __has_default_arithmetic_traits : public std::false_type {};


template <typename T1>
struct __has_default_arithmetic_traits<
        T1,
        void_t<decltype(DefaultArithmeticTraits<T1>::getTraits())>
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
struct IsTraits : public Impl::__is_traits<T> {};

template<typename T>
struct HasDefaultTraits : public Impl::__has_default_traits<T> {};

template<typename T>
struct HasDefaultIOTraits : public Impl::__has_default_io_traits<T> {};

template<typename T>
struct HasDefaultArithmeticTraits : public Impl::__has_default_arithmetic_traits<T> {};

#endif // CUSTOMTRAITS_H
