#ifndef MEMBERAPPROACH_H
#define MEMBERAPPROACH_H
#include "AccessType.h"
#include "GetSetAccess.h"



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

// Experimental function traits approach


template <class Class, class GetFunctor, class SetFunctor>
class MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<!std::is_bind_expression<GetFunctor>::value>>
: public Impl::GetAccess<GetFunctor>, public Impl::SetAccess<SetFunctor>{
public:
    using typeValue = std::decay_t<typename function_traits<GetFunctor>::return_type>;
    using typeClass = std::decay_t<typename function_traits<GetFunctor>::template argument<0>::type>;

public:

    MemberAccess(std::pair<GetFunctor, SetFunctor> getSet)
        : Impl::GetAccess<GetFunctor>(getSet.first), Impl::SetAccess<SetFunctor>(getSet.second)
    {}

    MemberAccess(const MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<!std::is_bind_expression<GetFunctor>::value>>&) = default;

    MemberAccess(MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<!std::is_bind_expression<GetFunctor>::value>>&&) = default;

};



template <class Class, class GetFunctor, class SetFunctor>
class MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<std::is_bind_expression<GetFunctor>::value> >
: public Impl::GetAccess<GetFunctor, Class>, public Impl::SetAccess<SetFunctor, Class>{
public:
    using typeValue = typename Impl::GetAccess<GetFunctor, Class>::_typeValue;
    using typeClass = Class;

public:

    MemberAccess(std::pair<GetFunctor, SetFunctor> getSet)
        : Impl::GetAccess<GetFunctor, Class>(getSet.first), Impl::SetAccess<SetFunctor, Class>(getSet.second)
    {}

    MemberAccess(const MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<std::is_bind_expression<GetFunctor>::value> >&) = default;

    MemberAccess(MemberAccess<Class, std::pair<GetFunctor, SetFunctor>, std::enable_if_t<std::is_bind_expression<GetFunctor>::value> >&&) = default;

};



#endif // MEMBERAPPROACH_H
