#ifndef MEMBERAPPROACH_H
#define MEMBERAPPROACH_H
#include <type_traits>
#include <utility>


/**
 * @brief The MemberApproach class
 * Generic abstract class providing the approach to
 * any attribute of a class using getValue and setValue
 */
template <typename Class, typename ValueType>
class MemberApproach{
public:
    virtual ValueType getValue(const Class*) const = 0;
    virtual void setValue(Class*, const ValueType&) const = 0;

    virtual ValueType getValue(const Class&) const = 0;
    virtual void setValue(Class&, const ValueType&) const = 0;
};


template<typename Class, typename ValueType, typename Ref>
class MemberReference{
    static_assert (std::is_same<Ref, ValueType Class::*>::value,
                   "The type MemberRef must be reference to member ValueType Class::* or pointer to getter and setter");
    MemberReference(Ref);
};


template <typename Class, typename ValueType>
class MemberReference<Class, ValueType, ValueType Class::*> : public MemberApproach<Class, ValueType>{

    using refType = ValueType Class::*;

    const refType ref;

public:

    MemberReference(refType referenceToMember) : ref(referenceToMember){
        //ref = referenceToMember;
    }

    MemberReference(const MemberReference<Class, ValueType, ValueType Class::*>&) = default;

    MemberReference(MemberReference<Class, ValueType, ValueType Class::*>&&) = default;

    virtual ValueType getValue(const Class* c) const override {
        return c->*ref;
    }

    virtual void setValue(Class* c, const ValueType& val) const override {
        c->*ref = val;
    }

    ValueType& getAttr(Class* c) const {
        return c->*ref;
    }

    virtual ValueType getValue(const Class& c) const override {
        return c.*ref;
    }

    virtual void setValue(Class& c, const ValueType& val) const override {
        c.*ref = val;
    }

    ValueType& getAttr(Class& c) const {
        return c.*ref;
    }
};





template <typename Class, typename ValueType>
class MemberReference<Class, ValueType, ValueType& (Class::*)()> : public MemberApproach<Class, ValueType>{

    using refType = ValueType& (Class::*)();

    refType const ref;

public:

    MemberReference(refType referenceToMember)
        : ref(referenceToMember)
    {
        //ref = referenceToMember;
    }

    virtual ValueType getValue(const Class* c) const override {
        return (c->*ref)();
    }

    virtual void setValue(Class* c, const ValueType& val) const override {
        (c->*ref)() = val;
    }
    virtual ValueType getValue(const Class& c) const override {
        return (c.*ref)();
    }

    virtual void setValue(Class& c, const ValueType& val) const override {
        (c.*ref)() = val;
    }

};




template <typename Class, typename ValueType>
class MemberReference<
        Class,
        ValueType,
        std::pair<ValueType (Class::*)(), void (Class::*)(const ValueType&)>
>
: public MemberApproach<Class, ValueType>{

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

     virtual ValueType getValue(const Class* c) const override {
         return (c->*refGet)();
     }

     virtual void setValue(Class* c, const ValueType& val) const override {
         (c->*refSet)(val);
     }

     virtual ValueType getValue(const Class& c) const override {
         return (c.*refGet)();
     }

     virtual void setValue(Class& c, const ValueType& val) const override {
         (c.*refSet)(val);
     }
};


template <typename Class, typename ValueType>
class MemberReference<
        Class,
        ValueType,
        std::pair<ValueType (Class::*)() const, void (Class::*)(const ValueType&)>
>
: public MemberApproach<Class, ValueType>{

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

     virtual ValueType getValue(const Class* c) const override {
         return (c->*refGet)();
     }

     virtual void setValue(Class* c, const ValueType& val) const override {
         (c->*refSet)(val);
     }

     virtual ValueType getValue(const Class& c) const override {
         return (c.*refGet)();
     }

     virtual void setValue(Class& c, const ValueType& val) const override {
         (c.*refSet)(val);
     }
};


template <typename Class, typename ValueType>
class MemberReference<
        Class,
        ValueType,
        std::pair<const ValueType& (Class::*)() const, void (Class::*)(const ValueType&)>
>
: public MemberApproach<Class, ValueType>{

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

     virtual ValueType getValue(const Class* c) const override {
         return (c->*refGet)();
     }

     virtual void setValue(Class* c, const ValueType& val) const override {
         (c->*refSet)(val);
     }

     virtual ValueType getValue(const Class& c) const override {
         return (c.*refGet)();
     }

     virtual void setValue(Class& c, const ValueType& val) const override {
         (c.*refSet)(val);
     }
};


template <typename Class, typename ValueType>
class MemberReference<
        Class,
        ValueType,
        std::pair<const ValueType& (Class::*)(), void (Class::*)(const ValueType&)>
>
: public MemberApproach<Class, ValueType>{

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

     virtual ValueType getValue(const Class* c) const override {
         return (c->*refGet)();
     }

     virtual void setValue(Class* c, const ValueType& val) const override {
         (c->*refSet)(val);
     }

     virtual ValueType getValue(const Class& c) const override {
         return (c.*refGet)();
     }

     virtual void setValue(Class& c, const ValueType& val) const override {
         (c.*refSet)(val);
     }
};



template <typename Class, typename ValueType>
class MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class&), const ValueType&(*)(const Class&)>> : public MemberApproach<Class, ValueType>{

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

    virtual ValueType getValue(const Class* c) const override {
        return get(*c);
    }


    virtual void setValue(Class* c, const ValueType& val) const override {
        set(*c) = val;
    }

    virtual ValueType getValue(const Class& c) const override {
        return get(c);
    }


    virtual void setValue(Class& c, const ValueType& val) const override {
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
class MemberReference<Class, ValueType, std::pair<ValueType&(*)(Class*), const ValueType&(*)(const Class*)>> : public MemberApproach<Class, ValueType>{

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

    virtual ValueType getValue(const Class* c) const override {
        return get(c);
    }


    virtual void setValue(Class* c, const ValueType& val) const override {
        set(c) = val;
    }

    virtual ValueType getValue(const Class& c) const override {
        return get(&c);
    }


    virtual void setValue(Class& c, const ValueType& val) const override {
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
