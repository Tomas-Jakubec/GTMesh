#ifndef INTERFACECOMPACT_H
#define INTERFACECOMPACT_H
#include <type_traits>
#include "../CustomTypeTraits.h"
namespace Interface {
/**
 * Generic template for indexable interface
 */
namespace Impl {
template<typename T>
constexpr bool IsSimpleSerializable_v = std::is_arithmetic<T>::value || std::is_enum<T>::value;


template<class T>
constexpr bool _has_constexpr_size(const T& val) { return val.size() == 0 || val.size() != 0; }
//template<class T, ::Impl::void_t<decltype (&T::getSize)>* = nullptr, int=T::getSize()>
//constexpr bool _has_constexpr_size(T) { return true; }
constexpr bool _has_constexpr_size(...) { return false; }

template<typename T>
class HasConstexprSize : public std::conditional<_has_constexpr_size(std::declval<T>()),
                                                 std::true_type,
                                                 std::false_type>
{};
} // namespace Impl

template<typename T, typename = void>
struct Compact : public std::false_type
{};

template<typename T>
struct Compact<T, std::enable_if_t<Impl::IsSimpleSerializable_v<T>>> : public std::true_type
{};

template<typename T>
struct Compact<T, std::enable_if_t<Impl::HasConstexprSize<T>::value>> : public std::true_type
{};

}
#endif // INTERFACECOMPACT_H
