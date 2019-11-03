#ifndef TRAITS_H
#define TRAITS_H
#include "MemberApproach.h"
#include "../debug/Debug.h"
#include <string>
#include <memory>
#include "Singleton.h"

template<typename Class, typename...Types>
class Traits {

    template <unsigned int Index>
    using type = std::tuple_element_t<Index,std::tuple<Types...>>;

    template<unsigned int Index = sizeof...(Types) - 1, typename Dummy = void>
    struct MemRefs: public MemRefs<Index - 1>{
        std::unique_ptr<MemberApproach<Class, type<Index>>> ref = nullptr;
        std::string name;
        MemRefs(){
            DBGVAR(this->name, Index);
        }
    };

    template<typename Dummy>
    struct MemRefs<0, Dummy>{
        std::unique_ptr<MemberApproach<Class, type<0>>> ref = nullptr;
        std::string name;
        MemRefs(){
            DBGVAR(this->name, 0);
        }
    };

    using refs = Singleton<MemRefs<sizeof... (Types) - 1, void>>;

    template<unsigned int Pos = 0, typename ref, typename...Refs>
    static void _makeReferences(const std::string& name, ref member,Refs... refsAndNames) {
        _makeReferences<Pos, ref>(name, member);
        _makeReferences<Pos+1>(refsAndNames...);
    }

    template<unsigned int Pos, typename ref>
    static void _makeReferences(const std::string& name, ref member) {
        refs::getInstance().MemRefs<Pos, void>::name = name;
        DBGVAR((refs::getInstance().MemRefs<Pos, void>::name));
        getReference<Pos>() = std::unique_ptr<MemberApproach<Class, type<Pos>>>(new typename MemberReference<Class, type<Pos>>::template SuperRef<ref>(member));
    }


    template<unsigned int Pos, typename ref>
    static void _makeReferences(const char* name, ref member) {
        refs::getInstance().MemRefs<Pos, void>::name = name;
        DBGVAR((refs::getInstance().MemRefs<Pos, void>::name));
        getReference<Pos>() = std::unique_ptr<MemberApproach<Class, type<Pos>>>(new typename MemberReference<Class, type<Pos>>::template SuperRef<ref>(member));
    }


public:



    static constexpr unsigned int size(){
        return sizeof... (Types);
    }


    template<unsigned int Index>
    static std::unique_ptr<MemberApproach<Class, type<Index>>>& getReference(){
        return refs::getInstance().MemRefs<Index, void>::ref;
    }

    template<unsigned int Index>
    static std::string& getName(){
        return refs::getInstance().MemRefs<Index, void>::name;
    }


    template<typename...Refs>
    Traits(Refs... refsAndNames){
        DBGMSG("making refs");
        makeReferences(refsAndNames...);
        DBGVAR(getName<0>(), getName<1>());
    }

    template<typename...Refs>
    static void makeReferences(Refs... refsAndNames) {
        _makeReferences<0>(refsAndNames...);
    }

};





template <typename Ref>
struct MemberReferenceType {
    static_assert (std::is_trivial_v<Ref>, "The Ref must be a type of member reference (MemberType Class::*) or reference getter (MemberType& (Class::*)()) into class \
or a pair of getter and setter (std::pair<MemberType (Class::*)(), void (Class::*)(const MemberType&))");
};

template <typename Class, typename MemberType>
struct MemberReferenceType <MemberType Class::*>{
    using type = MemberType;
};


template <typename Class, typename MemberType>
struct MemberReferenceType <MemberType& (Class::*)()>{
    using type = MemberType;
};


template <typename Class, typename MemberType>
struct MemberReferenceType <std::pair<MemberType (Class::*)(), void (Class::*)(const MemberType&)>>{
    using type = MemberType;
};


template<typename Class>
class Traits<Class>: public std::false_type{
    static_assert (true, "The Traits template must be specialized for given type and must contain Traits references using variadic Traits.");
};



#include "../../../Macros/MacroForEach.h"


#define MEMREF_TYPE_CUSTOM(name, memberRef) typename MemberReferenceType<decltype(memberRef)>::type
#define MAKE_CUSTOM_ATTRIBUTE_TRAIT(Class,...)\
template<>                              \
class Traits<Class>: public std::true_type{                 \
public:                             \
    using ttype = Traits<Class, FOR_EACH_2ARGS(MEMREF_TYPE_CUSTOM, __VA_ARGS__)>; \
    const static ttype tr;   \
    Traits(){DBGMSG("TRAITS");}\
}; \
const Traits<Class>::ttype Traits<Class>::tr(__VA_ARGS__);

#define NAME_AND_REF(Class, name, member) name, &Class::member

#define MAKE_NAMED_ATRIBUTE_TRAIT(Class, ...) MAKE_CUSTOM_ATTRIBUTE_TRAIT(Class, FOR_EACH_3ARGS_1STAT(NAME_AND_REF, Class, __VA_ARGS__))

#define NAME_ATT(attribute) #attribute, attribute
#define MAKE_ATRIBUTE_TRAIT(Class, ...) MAKE_NAMED_ATRIBUTE_TRAIT(Class, FOR_EACH(NAME_ATT, __VA_ARGS__))



#endif // TRAITS_H
