#ifndef MEMBERAPPROACH_H
#define MEMBERAPPROACH_H
#include <type_traits>
#include <utility>



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



template<typename Ref>
class MemberAccess{
    static_assert (!std::is_same<Ref, void>::value,
                   "The type Ref must be reference to member (ValueType Class::*),"
                   " member function or pointer to getter and setter");
    MemberAccess(Ref);
};


template <typename Class, typename ValueType>
class MemberAccess<ValueType Class::*> : public DirectAccess, public ConstGetAccess {
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

    MemberAccess(const MemberAccess<ValueType Class::*>&) = default;

    MemberAccess(MemberAccess<ValueType Class::*>&&) = default;

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

#endif // MEMBERAPPROACH_H
