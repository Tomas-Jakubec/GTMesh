#ifndef TRAITS_H
#define TRAITS_H
#include "MemberApproach.h"
#include <string>

namespace Detail {


}


template<typename Class, typename...Types>
class Traits {

    template <unsigned int Index>
    using type = std::tuple_element_t<Index,std::tuple<Types...>>;

    template<unsigned int Index = sizeof...(Types) - 1, typename Dummy = void>
    struct MemRefs: public MemRefs<Index - 1>{
        MemberApproach<Class, type<Index>>* ref = nullptr;
        std::string name;
    };

    template<typename Dummy>
    struct MemRefs<0, Dummy>{
        MemberApproach<Class, type<0>>* ref = nullptr;
        std::string name;
    };

public:
    static MemRefs<sizeof... (Types) - 1, void> refs;

    template<unsigned int Index>
    static MemberApproach<Class, type<Index>>*& getReference(){
        return refs.MemRefs<Index, void>::ref;
    }


    template<typename...Refs>
    static void makeReferences(Refs... refsAndNames) {
        _makeReferences<0>(refsAndNames...);
    }

    template<unsigned int Pos = 0, typename ref, typename...Refs>
    static void _makeReferences(std::string name, ref member,Refs... refsAndNames) {
        _makeReferences<Pos, ref>(name, member);
        _makeReferences<Pos+1, Refs...>(refsAndNames...);
    }

    template<unsigned int Pos, typename ref>
    static void _makeReferences(std::string name, ref member) {
        refs.MemRefs<Pos, void>::name = name;
        refs.MemRefs<Pos, void>::ref = new typename MemberReference<Class, type<Pos>>::Reference(member);
    }




};

template<typename Class, typename... Types>
typename Traits<Class, Types...>::template MemRefs<sizeof... (Types) - 1> Traits<Class,Types...>::refs;

#endif // TRAITS_H
