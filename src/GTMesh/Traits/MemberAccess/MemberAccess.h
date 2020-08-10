#ifndef MEMBERAPPROACH_H
#define MEMBERAPPROACH_H
#include <type_traits>
#include <utility>
#include <functional>


/**
 * @brief The DirectAccess struct determines that the
 * reference can provide direct approach to the member.
 */
struct DirectAccess{
    static constexpr std::true_type is_direct{};
};


/**
 * @brief The DirectAccess struct determines that the
 * reference can provide constant access to the member.
 */
struct ConstGetAccess{
    static constexpr std::true_type has_const_get{};
};


namespace Impl {
template <typename T, typename = void>
struct IsDirectAccess : public std::false_type {};

template <typename T>
struct IsDirectAccess
        <
        T,
        typename std::enable_if<T::is_direct>::type
        >
        : public std::true_type {};

template <typename T, typename = void>
struct HasConstGetAccess : public std::false_type {};

template <typename T>
struct HasConstGetAccess
        <
        T,
        typename std::enable_if<T::has_const_get>::type
        >
        : public std::true_type {};


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



template<typename Class, typename Ref, typename = void>
class MemberAccess{
    static_assert (!std::is_same<Ref, void>::value,
                   "The type Ref must be reference to member (ValueType Class::*),"
                   " member function or pointer to getter and setter");
public:
    MemberAccess(Ref);
};


template <typename Class, typename ValueType>
class MemberAccess<Class, ValueType Class::*, void> : public DirectAccess, public ConstGetAccess {
public:
    using typeValue = ValueType;
    using typeClass = Class;
    using refType = ValueType Class::*;

private:
    const refType ref;

public:

    MemberAccess(refType referenceToMember) : ref(referenceToMember){
        //ref = referenceToMember;
    }

    MemberAccess(const MemberAccess<Class, ValueType Class::*, void>&) = default;

    MemberAccess(MemberAccess<Class, ValueType Class::*, void>&&) = default;

    ValueType getValue(const Class* c) const {
        return c->*ref;
    }

    void setValue(Class* c, const ValueType& val) const {
        c->*ref = val;
    }

    ValueType& getAttr(Class* c) const {
        return c->*ref;
    }

    ValueType getValue(const Class& c) const {
        return c.*ref;
    }

    void setValue(Class& c, const ValueType& val) const {
        c.*ref = val;
    }

    ValueType& getAttr(Class& c) const {
        return c.*ref;
    }
};




/*
template <typename Class, typename ValueType>
class MemberAccess<ValueType& (Class::*)()> : public DirectAccess {
public:
    using typeValue = ValueType;
    using typeClass = Class;

    using refType = ValueType& (Class::*)();

private:
    refType const ref;

public:

    MemberAccess(refType referenceToMember)
        : ref(referenceToMember)
    {
        //ref = referenceToMember;
    }

    MemberAccess(const MemberAccess<refType>&) = default;

    MemberAccess(MemberAccess<refType>&&) = default;

    ValueType getValue(Class* c) const {
        return (c->*ref)();
    }

    ValueType& getAttr(Class* c) const {
        return (c->*ref)();
    }

    void setValue(Class* c, const ValueType& val) const {
        (c->*ref)() = val;
    }

    ValueType getValue(Class& c) const {
        return (c.*ref)();
    }

    void setValue(Class& c, const ValueType& val) const {
        (c.*ref)() = val;
    }

    ValueType& getAttr(Class& c) const {
        return (c->*ref)();
    }

};




template <typename Class, typename ValueType>
class MemberAccess<std::pair<ValueType (Class::*)(), void (Class::*)(const ValueType&)>>
{
public:
    using typeValue = ValueType;
    using typeClass = Class;
    using getterType = ValueType (Class::*)();
    using setterType = void (Class::*)(const ValueType&);
private:
    getterType const refGet;
    setterType const refSet;

public:

    MemberAccess(std::pair<getterType, setterType> getSet)
        :refGet(getSet.first), refSet(getSet.second){
        //refGet = getSet.first;
        //refSet = getSet.second;
    }

    MemberAccess(const MemberAccess<std::pair<getterType, setterType>>&) = default;

    MemberAccess(MemberAccess<std::pair<getterType, setterType>>&&) = default;

    ValueType getValue(Class* c) const {
        return (c->*refGet)();
    }

    void setValue(Class* c, const ValueType& val) const {
        (c->*refSet)(val);
    }

    ValueType getValue(Class& c) const {
        return (c.*refGet)();
    }

    void setValue(Class& c, const ValueType& val) const {
        (c.*refSet)(val);
    }
};


template <typename Class, typename ValueType>
class MemberAccess<std::pair<ValueType (Class::*)() const, void (Class::*)(const ValueType&)>>
    : public ConstGetAccess {
public:
    using typeValue = ValueType;
    using typeClass = Class;

    using getterType = ValueType (Class::*)() const;
    using setterType = void (Class::*)(const ValueType&);

private:
    getterType const refGet;
    setterType const refSet;

public:

    MemberAccess(std::pair<getterType, setterType> getSet)
        :refGet(getSet.first), refSet(getSet.second){
        //refGet = getSet.first;
        //refSet = getSet.second;
    }

    MemberAccess(const MemberAccess<std::pair<getterType, setterType>>&) = default;

    MemberAccess(MemberAccess<std::pair<getterType, setterType>>&&) = default;


    ValueType getValue(const Class* c) const {
        return (c->*refGet)();
    }

    void setValue(Class* c, const ValueType& val) const {
        (c->*refSet)(val);
    }

    ValueType getValue(const Class& c) const {
        return (c.*refGet)();
    }

    void setValue(Class& c, const ValueType& val) const {
        (c.*refSet)(val);
    }
};


template <typename Class, typename ValueType>
class MemberAccess<std::pair<const ValueType& (Class::*)() const, void (Class::*)(const ValueType&)>>
    : public ConstGetAccess {

public:
    using typeValue = ValueType;
    using typeClass = Class;

    using getterType = const ValueType& (Class::*)() const;
    using setterType = void (Class::*)(const ValueType&);
private:
    getterType const refGet;
    setterType const refSet;

public:

     MemberAccess(std::pair<getterType, setterType> getSet)
         :refGet(getSet.first), refSet(getSet.second){
         //refGet = getSet.first;
         //refSet = getSet.second;
     }
     MemberAccess(const MemberAccess<std::pair<getterType, setterType>>&) = default;

     MemberAccess(MemberAccess<std::pair<getterType, setterType>>&&) = default;

     ValueType getValue(const Class* c) const {
         return (c->*refGet)();
     }

     void setValue(Class* c, const ValueType& val) const {
         (c->*refSet)(val);
     }

     ValueType getValue(const Class& c) const {
         return (c.*refGet)();
     }

     void setValue(Class& c, const ValueType& val) const {
         (c.*refSet)(val);
     }
};


template <typename Class, typename ValueType>
class MemberAccess<std::pair<const ValueType& (Class::*)(), void (Class::*)(const ValueType&)>>
{
public:
    using typeValue = ValueType;
    using typeClass = Class;

    using getterType = const ValueType& (Class::*)();
    using setterType = void (Class::*)(const ValueType&);

private:
    getterType const refGet;
    setterType const refSet;

public:

    MemberAccess(std::pair<getterType, setterType> getSet)
        :refGet(getSet.first), refSet(getSet.second){
        //refGet = getSet.first;
        //refSet = getSet.second;
    }
    MemberAccess(const MemberAccess<std::pair<getterType, setterType>>&) = default;

    MemberAccess(MemberAccess<std::pair<getterType, setterType>>&&) = default;

     ValueType getValue(Class* c) const {
         return (c->*refGet)();
     }

     void setValue(Class* c, const ValueType& val) const {
         (c->*refSet)(val);
     }

     ValueType getValue(Class& c) const {
         return (c.*refGet)();
     }

     void setValue(Class& c, const ValueType& val) const {
         (c.*refSet)(val);
     }
};



template <typename Class, typename ValueType>
class MemberAccess<std::pair<const ValueType&(*)(const Class&), ValueType&(*)(Class&)>>
    : public DirectAccess, public ConstGetAccess {
public:
    using typeValue = ValueType;
    using typeClass = Class;

    using setterType = ValueType& (*)(Class&);


    using getterType = const ValueType&(*)(const Class&);
private:
    const setterType set;
    const getterType get;

public:

    MemberAccess(std::pair<getterType, setterType> getSet)
        : set(getSet.second), get(getSet.first)
    {}

    MemberAccess(const MemberAccess<std::pair<getterType, setterType>>&) = default;

    MemberAccess(MemberAccess<std::pair<getterType, setterType>>&&) = default;

    ValueType getValue(const Class* c) const {
        return get(*c);
    }


    void setValue(Class* c, const ValueType& val) const {
        set(*c) = val;
    }

    ValueType getValue(const Class& c) const {
        return get(c);
    }


    void setValue(Class& c, const ValueType& val) const {
        set(c) = val;
    }

    ValueType& getAttr(Class& c) const {
        return set(c);
    }

    ValueType& getAttr(Class* c) const {
        return set(*c);
    }
};


template <typename Class, typename ValueType>
class MemberAccess<std::pair<const ValueType&(*)(const Class*), ValueType&(*)(Class*)>>
: public DirectAccess, public ConstGetAccess{
public:
    using typeValue = ValueType;
    using typeClass = Class;

    using setterType = ValueType& (*)(Class*);


    using getterType = const ValueType&(*)(const Class*);
private:
    const setterType set;
    const getterType get;

public:

    MemberAccess(std::pair<getterType, setterType> getSet)
        : set(getSet.second), get(getSet.first)
    {}

    MemberAccess(const MemberAccess<std::pair<getterType, setterType>>&) = default;

    MemberAccess(MemberAccess<std::pair<getterType, setterType>>&&) = default;



    void setValue(Class* c, const ValueType& val) const {
        set(c) = val;
    }


    ValueType getValue(const Class* c) const {
        return get(c);
    }

    ValueType getValue(const Class& c) const {
        return get(&c);
    }


    void setValue(Class& c, const ValueType& val) const {
        set(&c) = val;
    }

    ValueType& getAttr(Class& c) const {
        return set(&c);
    }

    ValueType& getAttr(Class* c) const {
        return set(c);
    }
};
*/
// Experimental function traits approach
#include "../FunctionTraits.h"
#include <type_traits>
namespace Impl {


struct empty{};

template<class GetFunctor, typename Class = void>
struct _get: public ConstGetAccess{
#if __cplusplus <= 201702L // standard c++14 and older
    using _typeValue = typename std::result_of<GetFunctor(Class)>::type;
#else
    using _typeValue = typename std::invoke_result<GetFunctor,Class>::type;
#endif
    using _getterType = GetFunctor;
private:
    const _getterType get;
public:

    _get(_getterType& g): get(g){};

    _get(const _get< GetFunctor, Class >&) = default;

    _get(_get< GetFunctor, Class >&&) = default;


    auto getValue(const Class* c) const {
        return get(*c);
    }

    auto getValue(const Class& c) const {
        return get(c);
    }
};

template <class GetFunctor>
struct _get<GetFunctor, std::enable_if_t<!std::is_member_pointer<GetFunctor>::value && !std::is_bind_expression<GetFunctor>::value, void>>
: public std::conditional_t< std::is_const< std::remove_reference_t<typename function_traits< GetFunctor >::template argument< 0 >::type > >::value, ConstGetAccess, empty >
{

    using _getterType = GetFunctor;
private:
    const _getterType get;
    using _getClassRef = std::remove_reference_t<typename function_traits<GetFunctor>::template argument<0>::type>;
public:

    _get(_getterType& g): get(g){};

    _get(const _get<GetFunctor, std::enable_if_t<!std::is_member_pointer<GetFunctor>::value && !std::is_bind_expression<GetFunctor>::value, void>>&) = default;

    _get(_get<GetFunctor, std::enable_if_t<!std::is_member_pointer<GetFunctor>::value && !std::is_bind_expression<GetFunctor>::value, void>>&&) = default;


    auto getValue(_getClassRef* c) const {
        return get(*c);
    }

    auto getValue(_getClassRef& c) const {
        return get(c);
    }
};


template <class GetFunctor>
struct _get<GetFunctor, std::enable_if_t<std::is_member_pointer<GetFunctor>::value, void>>
: public std::conditional_t< std::is_const< std::remove_reference_t<typename function_traits< GetFunctor >::template argument< 0 >::type > >::value, ConstGetAccess, empty >
{
private:

    using _getterType = GetFunctor;
    const _getterType get;
    using _getClassRef = std::remove_reference_t<typename function_traits<GetFunctor>::template argument<0>::type>;
public:

    _get(_getterType& g): get(g){};

    _get(const _get<GetFunctor, std::enable_if_t<std::is_member_pointer<GetFunctor>::value, void>>&) = default;

    _get(_get<GetFunctor, std::enable_if_t<std::is_member_pointer<GetFunctor>::value, void>>&&) = default;


    auto getValue(_getClassRef* c) const {
        return (c->*get)();
    }

    /**
     * @brief Returns a value of a member of the object c.
     */
    auto getValue(_getClassRef& c) const {
        return (c.*get)();
    }
};


template <class SetFunctor, typename Spec = void>
struct _set{
    using _setterType = SetFunctor;
private:
    const _setterType set;
public:
    _set(const _setterType& s): set(s){}


    _set(const _set<SetFunctor, Spec>&) = default;

    _set(_set<SetFunctor, Spec>&&) = default;

    template<typename typeClass, typename typeValue>
    void setValue(typeClass* c, const typeValue& val) const {
        set(c, val);
    }

    template<typename typeClass, typename typeValue>
    void setValue(typeClass& c, const typeValue& val) const {
        set(c, val);
    }
};


template<typename SetFunctor>
struct _set<SetFunctor, std::enable_if_t<!std::is_member_pointer<SetFunctor>::value && function_traits<SetFunctor>::arity == 2, void>>{

    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::template argument<1>::type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;
private:
    const _setterType set;
public:
    _set(const _setterType& s): set(s){}


    void setValue(_typeClass* c, const _typeValue& val) const {
        set(c) = val;
    }


    void setValue(_typeClass& c, const _typeValue& val) const {
        set(c) = val;
    }

    auto getAttr(_typeClass* c) const {
        return set(c);
    }


    auto getAttr(_typeClass& c) const {
        return set(c);
    }
};


template<typename SetFunctor>
struct _set<SetFunctor, std::enable_if_t<!std::is_member_pointer<SetFunctor>::value && function_traits<SetFunctor>::arity == 1, void>>{

    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::return_type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;
private:
    const _setterType set;
public:
    _set(const _setterType& s): set(s){}



    void setValue(_typeClass* c, const _typeValue& val) const {
        set(c) = val;
    }


    void setValue(_typeClass& c, const _typeValue& val) const {
        set(c) = val;
    }
};

template<typename SetFunctor>
struct _set<SetFunctor, std::enable_if_t<std::is_member_pointer<SetFunctor>::value && function_traits<SetFunctor>::arity == 2, void>>{
private:
    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::template argument<1>::type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;
private:
    const _setterType set;
public:
    _set(const _setterType& s): set(s){}

    void setValue(_typeClass* c, const _typeValue& val) const {
        (c->*set)(val);
    }


    void setValue(_typeClass& c, const _typeValue& val) const {
        (c.*set)(val);
    }
};

template<typename SetFunctor>
struct _set<SetFunctor, std::enable_if_t<std::is_member_pointer<SetFunctor>::value && function_traits<SetFunctor>::arity == 1, void>>
: public DirectAccess {
private:
    using _typeValue = std::decay_t<typename function_traits<SetFunctor>::return_type>;
    using _typeClass = std::decay_t<typename function_traits<SetFunctor>::template argument<0>::type>;

    using _setterType = SetFunctor;
private:
    const _setterType set;
public:
    _set(const _setterType& s): set(s){}



    void setValue(_typeClass* c, const _typeValue& val) const {
        (c->*set)() = val;
    }


    void setValue(_typeClass& c, const _typeValue& val) const {
        (c.*set)() = val;
    }

    auto getAttr(_typeClass* c) const {
        return (c->*set)();
    }


    auto getAttr(_typeClass& c) const {
        return (c.*set)();
    }
};
}

template <class Class, class GetFunctor, class SetFunctor>
class MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<!std::is_bind_expression<GetFunctor>::value>>
: public Impl::_get<GetFunctor>, public Impl::_set<SetFunctor>{
public:
    using typeValue = std::decay_t<typename function_traits<GetFunctor>::return_type>;
    using typeClass = std::decay_t<typename function_traits<GetFunctor>::template argument<0>::type>;

public:

    MemberAccess(std::pair<GetFunctor, SetFunctor> getSet)
        : Impl::_get<GetFunctor>(getSet.first), Impl::_set<SetFunctor>(getSet.second)
    {}

    MemberAccess(const MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<!std::is_bind_expression<GetFunctor>::value>>&) = default;

    MemberAccess(MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<!std::is_bind_expression<GetFunctor>::value>>&&) = default;

};



template <class Class, class GetFunctor, class SetFunctor>
class MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<std::is_bind_expression<GetFunctor>::value> >
: public Impl::_get<GetFunctor, Class>, public Impl::_set<SetFunctor, Class>{
public:
    using typeValue = typename Impl::_get<GetFunctor, Class>::_typeValue;
    using typeClass = Class;

public:

    MemberAccess(std::pair<GetFunctor, SetFunctor> getSet)
        : Impl::_get<GetFunctor, Class>(getSet.first), Impl::_set<SetFunctor, Class>(getSet.second)
    {}

    MemberAccess(const MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<std::is_bind_expression<GetFunctor>::value> >&) = default;

    MemberAccess(MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<std::is_bind_expression<GetFunctor>::value> >&&) = default;

};



#endif // MEMBERAPPROACH_H
