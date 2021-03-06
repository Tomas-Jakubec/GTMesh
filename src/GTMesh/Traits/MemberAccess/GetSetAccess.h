#ifndef GETSETACCESS_H
#define GETSETACCESS_H
#include "../FunctionTraits.h"
#include "AccessType.h"
#include <functional>
#include <type_traits>
namespace Impl {

struct empty
{};

template<class GetFunctor, typename Class = void>
struct GetAccess : public ConstGetAccess
{
#if __cplusplus <= 201702L // standard c++14 and older
    using _typeValue = typename std::result_of<GetFunctor(Class)>::type;
#else
    using _typeValue = typename std::invoke_result<GetFunctor, Class>::type;
#endif
    using _getterType = GetFunctor;

private:
    const _getterType get;

public:
    constexpr GetAccess(const _getterType &g) : get(g){};

    constexpr GetAccess(const GetAccess<GetFunctor, Class> &) = default;

    constexpr GetAccess(GetAccess<GetFunctor, Class> &&) = default;

    inline auto getValue(const Class *c) const { return get(*c); }

    inline auto getValue(const Class &c) const { return get(c); }
};

template<class GetFunctor>
struct GetAccess<GetFunctor,
                 std::enable_if_t<!std::is_member_pointer<GetFunctor>::value
                                      && !std::is_bind_expression<GetFunctor>::value,
                                  void>>
    : public std::conditional_t<std::is_const<std::remove_reference_t<typename function_traits<
                                    GetFunctor>::template argument<0>::type>>::value,
                                ConstGetAccess,
                                empty>
{
    using _getterType = GetFunctor;

private:
    const _getterType get;
    using _getClassRef
        = std::remove_reference_t<typename function_traits<GetFunctor>::template argument<0>::type>;

public:
    constexpr GetAccess(const _getterType &g) : get(g){};

    constexpr GetAccess(
        const GetAccess<GetFunctor,
                        std::enable_if_t<!std::is_member_pointer<GetFunctor>::value
                                             && !std::is_bind_expression<GetFunctor>::value,
                                         void>> &)
        = default;

    constexpr GetAccess(GetAccess<GetFunctor,
                                  std::enable_if_t<!std::is_member_pointer<GetFunctor>::value
                                                       && !std::is_bind_expression<GetFunctor>::value,
                                                   void>> &&)
        = default;

    inline auto getValue(_getClassRef *c) const { return get(*c); }

    inline auto getValue(_getClassRef &c) const { return get(c); }
};

template<class GetFunctor>
struct GetAccess<GetFunctor, std::enable_if_t<std::is_member_pointer<GetFunctor>::value, void>>
    : public std::conditional_t<std::is_const<std::remove_reference_t<typename function_traits<
                                    GetFunctor>::template argument<0>::type>>::value,
                                ConstGetAccess,
                                empty>
{
private:
    using _getterType = GetFunctor;
    const _getterType get;
    using _getClassRef
        = std::remove_reference_t<typename function_traits<GetFunctor>::template argument<0>::type>;

public:
    constexpr GetAccess(const _getterType &g) : get(g){};

    constexpr GetAccess(
        const GetAccess<GetFunctor,
                        std::enable_if_t<std::is_member_pointer<GetFunctor>::value, void>> &)
        = default;

    constexpr GetAccess(
        GetAccess<GetFunctor, std::enable_if_t<std::is_member_pointer<GetFunctor>::value, void>> &&)
        = default;

    inline auto getValue(_getClassRef *c) const { return (c->*get)(); }

    /**
     * @brief Returns a value of a member of the object c.
     */
    inline auto getValue(_getClassRef &c) const { return (c.*get)(); }
};

template<class SetFunctor, typename Spec = void>
struct SetAccess
{
    using _setterType = SetFunctor;

private:
    const _setterType set;

public:
    constexpr SetAccess(const _setterType &s) : set(s) {}

    constexpr SetAccess(const SetAccess<SetFunctor, Spec> &) = default;

    constexpr SetAccess(SetAccess<SetFunctor, Spec> &&) = default;

    template<typename typeClass, typename typeValue>
    inline void setValue(typeClass *c, const typeValue &val) const
    {
        set(c, val);
    }

    template<typename typeClass, typename typeValue>
    inline void setValue(typeClass &c, const typeValue &val) const
    {
        set(c, val);
    }
};

template<typename SetFunctor>
struct SetAccess<SetFunctor,
                 std::enable_if_t<!std::is_member_pointer<SetFunctor>::value
                                      && function_traits<SetFunctor>::arity == 2,
                                  void>>
{
    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::template argument<1>::type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;

private:
    const _setterType set;

public:
    constexpr SetAccess(const _setterType &s) : set(s) {}

    inline void setValue(_typeClass *c, const _typeValue &val) const { set(c, val); }

    inline void setValue(_typeClass &c, const _typeValue &val) const { set(c, val); }
};

template<typename SetFunctor>
struct SetAccess<SetFunctor,
                 std::enable_if_t<!std::is_member_pointer<SetFunctor>::value
                                      && function_traits<SetFunctor>::arity == 1,
                                  void>> : public DirectAccess
{
    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::return_type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;

private:
    const _setterType set;

public:
    constexpr SetAccess(const _setterType &s) : set(s) {}

    inline void setValue(_typeClass *c, const _typeValue &val) const { set(c) = val; }

    inline void setValue(_typeClass &c, const _typeValue &val) const { set(c) = val; }

    inline auto &getAttr(_typeClass *c) const { return set(c); }

    inline auto &getAttr(_typeClass &c) const { return set(c); }
};

template<typename SetFunctor>
struct SetAccess<SetFunctor,
                 std::enable_if_t<std::is_member_pointer<SetFunctor>::value
                                      && function_traits<SetFunctor>::arity == 2,
                                  void>>
{
private:
    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::template argument<1>::type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;

private:
    const _setterType set;

public:
    constexpr SetAccess(const _setterType &s) : set(s) {}

    inline void setValue(_typeClass *c, const _typeValue &val) const { (c->*set)(val); }

    inline void setValue(_typeClass &c, const _typeValue &val) const { (c.*set)(val); }
};

template<typename SetFunctor>
struct SetAccess<SetFunctor,
                 std::enable_if_t<std::is_member_pointer<SetFunctor>::value &&
                                  function_traits<SetFunctor>::arity == 1,
                                  void>> : public DirectAccess
{
private:
    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::return_type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;

private:
    const _setterType set;

public:
    constexpr SetAccess(const _setterType &s) : set(s) {}

    inline void setValue(_typeClass *c, const _typeValue &val) const { (c->*set)() = val; }

    inline void setValue(_typeClass &c, const _typeValue &val) const { (c.*set)() = val; }

    inline auto &getAttr(_typeClass *c) const { return (c->*set)(); }

    inline auto &getAttr(_typeClass &c) const { return (c.*set)(); }
};
} // namespace Impl
#endif // GETSETACCESS_H
