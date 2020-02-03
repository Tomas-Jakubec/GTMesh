#ifndef MEMBERAPPROACH_H
#define MEMBERAPPROACH_H
#include <type_traits>
#include <utility>



/**
 * @brief The DirectReference struct determines that the
 * reference can provide direct approach to the member.
 */
struct DirectReference{
    static constexpr std::true_type is_direct{};
};





namespace Impl {
template <typename T, typename = void>
struct IsDirectReference : public std::false_type {};

template <typename T>
struct IsDirectReference
        <
        T,
        typename std::enable_if<T::is_direct>::type
        >
        : public std::true_type {};

} // Impl



/**
 * @brief The IsDirectReference struct inherits
 * @ref std::true_type if the class MemberReference provides direct
 * approach to the member using function getAttr.
 */
template <typename T>
struct IsDirectReference : public Impl::IsDirectReference<T> {};



template<typename Class, typename ValueType, typename Ref>
class MemberReference{
    static_assert (std::is_same<Ref, ValueType Class::*>::value,
                   "The type MemberRef must be reference to member ValueType Class::* or pointer to getter and setter");
    MemberReference(Ref);
};


template <typename Class, typename ValueType>
class MemberReference<Class, ValueType, ValueType Class::*> : public DirectReference {

    using refType = ValueType Class::*;

    const refType ref;

public:

    MemberReference(refType referenceToMember) : ref(referenceToMember){
        //ref = referenceToMember;
    }

    MemberReference(const MemberReference<Class, ValueType, ValueType Class::*>&) = default;

    MemberReference(MemberReference<Class, ValueType, ValueType Class::*>&&) = default;

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
class MemberReference<Class, ValueType, ValueType& (Class::*)()> : public DirectReference {

    using refType = ValueType& (Class::*)();

    refType const ref;

public:

    MemberReference(refType referenceToMember)
        : ref(referenceToMember)
    {
        //ref = referenceToMember;
    }

    ValueType getValue(const Class* c) const {
        return (c->*ref)();
    }

    ValueType& getAttr(Class* c) const {
        return (c->*ref)();
    }

    void setValue(Class* c, const ValueType& val) const {
        (c->*ref)() = val;
    }
    ValueType getValue(const Class& c) const {
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
class MemberReference<
        Class,
        ValueType,
        std::pair<ValueType (Class::*)(), void (Class::*)(const ValueType&)>
>
{

    using getter = ValueType (Class::*)();
    getter const refGet;
    using setter = void (Class::*)(const ValueType&);
    setter const refSet;

public:

     MemberReference(std::pair<getter, setter> getSet)
         :refGet(getSet.first), refSet(getSet.second){
         //refGet = getSet.first;
         //refSet = getSet.second;
     }

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
class MemberReference<
        Class,
        ValueType,
        std::pair<ValueType (Class::*)() const, void (Class::*)(const ValueType&)>
>
{

    using getter = ValueType (Class::*)() const;
    getter const refGet;
    using setter = void (Class::*)(const ValueType&);
    setter const refSet;

public:

     MemberReference(std::pair<getter, setter> getSet)
         :refGet(getSet.first), refSet(getSet.second){
         //refGet = getSet.first;
         //refSet = getSet.second;
     }

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
class MemberReference<
        Class,
        ValueType,
        std::pair<const ValueType& (Class::*)() const, void (Class::*)(const ValueType&)>
>
{

    using getter = const ValueType& (Class::*)() const;
    getter const refGet;
    using setter = void (Class::*)(const ValueType&);
    setter const refSet;

public:

     MemberReference(std::pair<getter, setter> getSet)
         :refGet(getSet.first), refSet(getSet.second){
         //refGet = getSet.first;
         //refSet = getSet.second;
     }

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
class MemberReference<
        Class,
        ValueType,
        std::pair<const ValueType& (Class::*)(), void (Class::*)(const ValueType&)>
>
{

    using getter = const ValueType& (Class::*)();
    getter const refGet;
    using setter = void (Class::*)(const ValueType&);
    setter const refSet;

public:

     MemberReference(std::pair<getter, setter> getSet)
         :refGet(getSet.first), refSet(getSet.second){
         //refGet = getSet.first;
         //refSet = getSet.second;
     }

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
class MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class&), const ValueType&(*)(const Class&)>>
: public DirectReference {

    using setter = ValueType& (*)(Class&);

    const setter set;

    using getter = const ValueType&(*)(const Class&);

    const getter get;

public:

    MemberReference(std::pair<ValueType&(*)(Class&), const ValueType&(*)(const Class&)> referenceToMember)
        : set(referenceToMember.first), get(referenceToMember.second)
    {
        //ref = referenceToMember;
    }

    MemberReference(const MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class&), const ValueType&(*)(const Class&)>>&) = default;

    MemberReference(MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class&), const ValueType&(*)(const Class&)>>&&) = default;

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
class MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class*), const ValueType&(*)(const Class*)>>
: public DirectReference {

    using setter = ValueType& (*)(Class*);

    const setter set;

    using getter = const ValueType&(*)(const Class*);

    const getter get;

public:

    MemberReference(std::pair<ValueType&(*)(Class*), const ValueType&(*)(const Class*)> referenceToMember)
        : set(referenceToMember.first), get(referenceToMember.second)
    {
        //ref = referenceToMember;
    }

    MemberReference(const MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class*), const ValueType&(*)(const Class*)>>&) = default;

    MemberReference(MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class*), const ValueType&(*)(const Class*)>>&&) = default;

    ValueType getValue(const Class* c) const {
        return get(c);
    }


    void setValue(Class* c, const ValueType& val) const {
        set(c) = val;
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





template <typename Ref>
struct MemberReferenceType {
    static_assert (std::is_trivial<Ref>::value, "The Ref must be a type of member reference (MemberType Class::*) or reference getter (MemberType& (Class::*)()) into class \
or a pair of getter and setter (std::pair<MemberType (Class::*)(), void (Class::*)(const MemberType&))");
};

template <typename Class, typename MemberType>
struct MemberReferenceType <MemberType Class::*>{
    using type = MemberType;
    using typeClass = Class;
};


template <typename Class, typename MemberType>
struct MemberReferenceType <MemberType& (Class::*)()>{
    using type = MemberType;
    using typeClass = Class;
};


template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<MemberType (Class::*)(), void (Class::*)(const MemberType&)>>{
    using type = MemberType;
    using typeClass = Class;
};

template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<MemberType (Class::*)()const, void (Class::*)(const MemberType&)>>{
    using type = MemberType;
    using typeClass = Class;
};

template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<const MemberType& (Class::*)()const, void (Class::*)(const MemberType&)>>{
    using type = MemberType;
    using typeClass = Class;
};

template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<const MemberType& (Class::*)(), void (Class::*)(const MemberType&)>>{
    using type = MemberType;
    using typeClass = Class;
};

template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<MemberType&(*)(Class&), const MemberType&(*)(const Class&)>>{
    using type = MemberType;
    using typeClass = Class;
};

template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<MemberType&(*)(Class*), const MemberType&(*)(const Class*)>>{
    using type = MemberType;
    using typeClass = Class;
};

#endif // MEMBERAPPROACH_H
