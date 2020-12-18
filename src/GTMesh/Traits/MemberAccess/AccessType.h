#ifndef ACCESSTYPE_H
#define ACCESSTYPE_H
#include <type_traits>
/**
 * @brief The DirectAccess struct determines that the
 * reference can provide direct approach to the member.
 */
struct DirectAccess{
    static constexpr bool is_direct = true;
};


/**
 * @brief The DirectAccess struct determines that the
 * reference can provide constant access to the member.
 */
struct ConstGetAccess{
    static constexpr bool has_const_get = true;
};


namespace Impl {
template <typename T, typename = void>
struct IsDirectAccess : public std::false_type {};

template<typename T>
struct IsDirectAccess<T, typename std::enable_if<T::is_direct>::type> : public std::true_type
{};

template <typename T, typename = void>
struct HasConstGetAccess : public std::false_type {};

template<typename T>
struct HasConstGetAccess<T, typename std::enable_if<T::has_const_get>::type> : public std::true_type
{};

} // Impl



/**
 * @brief The IsDirectReference struct inherits
 * @ref std::true_type if the class MemberReference provides direct
 * approach to the member using function getAttr.
 */
template <typename T>
struct IsDirectAccess : public Impl::IsDirectAccess<T> {};

template <typename T>
struct HasConstGetAccess : public Impl::HasConstGetAccess<T> {};

#endif // ACCESSTYPE_H
