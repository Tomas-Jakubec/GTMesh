#ifndef TRAITS_H
#define TRAITS_H
#include "MemberApproach.h"
#include <string>
#include <memory>
#include "Singleton.h"

template<typename Class, typename...Types>
class Traits {
public:
    template <unsigned int Index>
    using type = std::tuple_element_t<Index,std::tuple<Types...>>;
private:
    template<unsigned int Index = sizeof...(Types) - 1, typename Dummy = void>
    struct MemRefs: public MemRefs<Index - 1> {
        friend class Singleton<MemRefs<1, void>>;
        std::unique_ptr<MemberApproach<Class, type<Index>>> ref = nullptr;
        std::string name;

        MemRefs(){}
    };

    template<typename Dummy>
    struct MemRefs<0, Dummy>{
        friend class Singleton<MemRefs<0, void>>;
        std::unique_ptr<MemberApproach<Class, type<0>>> ref = nullptr;
        std::string name;

        MemRefs(){}
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
        refs::getInstance().MemRefs<Pos, void>::ref = std::unique_ptr<MemberApproach<Class, type<Pos>>>(new MemberReference<Class, type<Pos>, decltype(member)>(member));
    }


    template<unsigned int Pos, typename ref>
    static void _makeReferences(const char* name, ref member) {
        refs::getInstance().MemRefs<Pos, void>::name = name;
        refs::getInstance().MemRefs<Pos, void>::ref = std::unique_ptr<MemberApproach<Class, type<Pos>>>(new MemberReference<Class, type<Pos>, decltype(member)>(member));
    }


public:



    static constexpr unsigned int size(){
        return sizeof... (Types);
    }


    template<unsigned int Index>
    static const std::unique_ptr<MemberApproach<Class, type<Index>>>& getReference(){
        return refs::getInstance().MemRefs<Index, void>::ref;
    }

    template<unsigned int Index>
    static type<Index> getValue(const Class* c){
        return getReference<Index>()->getValue(c);
    }

    template<unsigned int Index>
    static type<Index> getValue(const Class& c){
        return getReference<Index>()->getValue(c);
    }

    template<unsigned int Index>
    static void setValue(Class* c, const type<Index>& val){
        getReference<Index>()->setValue(c, val);
    }

    template<unsigned int Index>
    static void setValue(Class& c, const type<Index>& val){
        getReference<Index>()->setValue(c, val);
    }


    template<unsigned int Index>
    static const std::string& getName(){
        return refs::getInstance().MemRefs<Index, void>::name;
    }


    template<typename...Refs>
    Traits(Refs... refsAndNames){
        makeReferences(refsAndNames...);
    }

    template<typename...Refs>
    static void makeReferences(Refs... refsAndNames) {
        _makeReferences<0>(refsAndNames...);
    }

};







template<typename Class>
class Traits<Class>: public std::false_type{
    static_assert (true, "The Traits template must be specialized for given type and must contain Traits references using variadic Traits.");
public:
    static constexpr std::false_type is_specialized{};
};

#include "../../Macros/MacroForEach.h"


#define MEMREF_TYPE_CUSTOM(name, memberRef) typename MemberReferenceType<decltype(memberRef)>::type
#define MAKE_CUSTOM_ATTRIBUTE_TRAIT(Class,...)\
template<>                              \
class Traits<Class>{                 \
public:                             \
    static constexpr std::true_type is_specialized{}; \
    using ttype = Traits<Class, FOR_EACH_2ARGS(MEMREF_TYPE_CUSTOM, __VA_ARGS__)>; \
    const static ttype tr;   \
    Traits(){DBGMSG("TRAITS");}\
}; \
const Traits<Class>::ttype Traits<Class>::tr(__VA_ARGS__); \

#define NAME_AND_REF(Class, name, member) name, &Class::member

#define MAKE_NAMED_ATTRIBUTE_TRAIT(Class, ...) MAKE_CUSTOM_ATTRIBUTE_TRAIT(Class, FOR_EACH_3ARGS_1STAT(NAME_AND_REF, Class, __VA_ARGS__))

#define NAME_ATT(attribute) #attribute, attribute
#define MAKE_ATTRIBUTE_TRAIT(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT(Class, FOR_EACH(NAME_ATT, __VA_ARGS__))



#endif // TRAITS_H
