#ifndef MEMBERAPPROACH_H
#define MEMBERAPPROACH_H
#include <type_traits>
#include <utility>

template <typename Class, typename ValueType>
class MemberApproach{
public:
    virtual ValueType getValue(Class*) = 0;
    virtual void setValue(Class*, const ValueType&) = 0;
};





template<typename Class, typename ValueType>
class MemberReference{

public:
    template<typename Ref, typename VOID = void>
    class Reference{
        static_assert (std::is_same<Ref, ValueType Class::*>::value, "The type MemberRef must be reference to member ValueType Class::* or pointer to getter and setter");

    public:
        Reference(Ref){}
    };


    template<typename Ref>
    class SuperRef: public Reference<Ref, void>{
    public:
        SuperRef(Ref ref): Reference<Ref, void> (ref){}
    };


   template<typename Ref>
   class Reference<
           Ref,
           // enable if MemberRef is pointer to member
           typename std::enable_if<
               std::is_same<Ref, ValueType Class::*>::value
           >::type
           >
   : public MemberApproach<Class, ValueType>{

       Ref ref;

   public:

       Reference(Ref referenceToMember){
           ref = referenceToMember;
           DBGVAR((std::is_same<Ref, ValueType Class::*>::value));
       }

       virtual ValueType getValue(Class* c) override {
           return c->*ref;
       }

       virtual void setValue(Class* c, const ValueType& val) override {
           c->*ref = val;
       }
   };







    template<typename Ref>
    class Reference<
            Ref,
            // enable if MemberRef is pointer to member
            typename std::enable_if<
                std::is_same<Ref, ValueType& (Class::*)()>::value
            >::type
            >
    : public MemberApproach<Class, ValueType>{

        Ref ref;

    public:

        Reference(Ref referenceToMember){
            ref = referenceToMember;
            DBGVAR((std::is_same<Ref, ValueType Class::*>::value));
        }

        virtual ValueType getValue(Class* c) override {
            return (c->*ref)();
        }

        virtual void setValue(Class* c, const ValueType& val) override {
            (c->*ref)() = val;
        }
    };




    template<typename MemberRefGet, typename MemberRefSet>
    class Reference<
            std::pair<MemberRefGet, MemberRefSet>,
            // enable if MemberRef is pointer to member
            typename std::enable_if<
                std::is_same<MemberRefGet, ValueType (Class::*)()>::value &&
                std::is_same<MemberRefSet, void (Class::*)(const ValueType&)>::value
            >::type
            >
    : public MemberApproach<Class, ValueType>{

        MemberRefGet refGet;
        MemberRefSet refSet;

    public:

        Reference(std::pair<MemberRefGet, MemberRefSet> getSet){
            refGet = getSet.first;
            refSet = getSet.second;
        }

        virtual ValueType getValue(Class* c) override {
            return (c->*refGet)();
        }

        virtual void setValue(Class* c, const ValueType& val) override {
            (c->*refSet)(val);
        }
    };

};
#endif // MEMBERAPPROACH_H
