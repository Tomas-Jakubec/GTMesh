#ifndef TRAITS_H
#define TRAITS_H
#include "MemberApproach.h"
#include <string>
#include <memory>


template<typename Class, typename...Types>
class Traits {

    template <unsigned int Index>
    using type = std::tuple_element_t<Index,std::tuple<Types...>>;

    template<unsigned int Index = sizeof...(Types) - 1, typename Dummy = void>
    struct MemRefs: public MemRefs<Index - 1>{
        std::unique_ptr<MemberApproach<Class, type<Index>>> ref = nullptr;
        std::string name;
    };

    template<typename Dummy>
    struct MemRefs<0, Dummy>{
        std::unique_ptr<MemberApproach<Class, type<0>>> ref = nullptr;
        std::string name;
    };



    template<unsigned int Pos = 0, typename ref, typename...Refs>
    static void _makeReferences(std::string name, ref member,Refs... refsAndNames) {
        _makeReferences<Pos, ref>(name, member);
        _makeReferences<Pos+1>(refsAndNames...);
    }

    template<unsigned int Pos, typename ref>
    static void _makeReferences(std::string name, ref member) {
        refs.MemRefs<Pos, void>::name = name;
        getReference<Pos>() = std::unique_ptr<MemberApproach<Class, type<Pos>>>(new typename MemberReference<Class, type<Pos>>::template SuperRef<ref>(member));
    }


public:
    static MemRefs<sizeof... (Types) - 1, void> refs;


    static constexpr unsigned int size(){
        return sizeof... (Types);
    }


    template<unsigned int Index>
    static std::unique_ptr<MemberApproach<Class, type<Index>>>& getReference(){
        return refs.MemRefs<Index, void>::ref;
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

template<typename Class, typename... Types>
typename Traits<Class, Types...>::template MemRefs<sizeof... (Types) - 1> Traits<Class,Types...>::refs;




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
class Traits<Class>{
    static_assert (true, "The Traits template must be specialized for given type and must contain Traits references using variadic Traits.");
};




#define _MAKE_ATRIBUTE_TRAITS(Class, ...)    \
template<>                              \
class Traits<Class>{                 \
public:                                 \
    static Traits<Class, double, Vector<3,double>> tr;   \
};                                      \
Traits<Class, double, Vector<3,double>> Traits<Class>::   \
tr(__VA_ARGS__); \



#define MAKE_ATRIBUTE_TRAITS(Class, ...)    \
template<>                              \
class Traits<Class>{                 \
public:                                 \
    static Traits<Class, double, Vector<3,double>> tr;   \
};                                      \
Traits<Class, double, Vector<3,double>> Traits<Class>::   \
tr(__VA_ARGS__); \


#endif // TRAITS_H
