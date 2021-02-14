#ifndef INTERFACECOMPACT_H
#define INTERFACECOMPACT_H
#include <type_traits>
#include "../CustomTypeTraits.h"
namespace Interface {
/**
 * Generic template for compact interface
 */
namespace Impl {

template<typename T>
struct IsSimpleSerializable
    : public std::conditional<std::is_arithmetic<T>::value || std::is_enum<T>::value,
                              std::true_type,
                              std::false_type>::type
{};

template<typename T>
constexpr bool IsSimpleSerializable_v = IsSimpleSerializable<T>::value;

} // namespace Impl

template<typename T, typename = void>
struct Compact : public std::false_type
{};

template<typename T>
struct Compact<T, std::enable_if_t<Impl::IsSimpleSerializable_v<T>>> : public std::true_type
{};

}
#endif // INTERFACECOMPACT_H
