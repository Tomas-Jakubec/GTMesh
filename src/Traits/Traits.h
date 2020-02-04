#ifndef TRAITS_H
#define TRAITS_H
#include "MemberApproach/MemberApproach.h"
#include <string>
#include <memory>
#include <functional>

template<typename Class, typename...RefTypes>
class Traits {
public:
    template<unsigned int Index>
    using refType = typename std::tuple_element<Index,std::tuple<RefTypes...>>::type;

    template <unsigned int Index>
    using type = typename MemberReferenceType<refType<Index>>::type;

private:
    template<unsigned int Index>
    using memRefType = MemberReference<Class, type<Index>, refType<Index>>;

    template<unsigned int Index = 0, typename Dummy = void>
    struct MemRefs: public MemRefs<Index + 1> {

        const MemberReference<Class, type<Index>, refType<Index>> ref;
        const char* name;

        template <typename ... REST>
        MemRefs(const char* n, refType<Index> r, REST... rest) : MemRefs<Index + 1> (rest...), ref(r), name(n){}
    };

    template<typename Dummy>
    struct MemRefs<sizeof...(RefTypes) - 1, Dummy>{
        const MemberReference<Class, type<sizeof...(RefTypes) - 1>, refType<sizeof...(RefTypes) - 1>> ref;
        const char* name;

        MemRefs(const char* n, refType<sizeof...(RefTypes) - 1> r) : ref(r), name(n){}
    };

    const MemRefs<0, void> refs;

public:



    static constexpr unsigned int size(){
        return sizeof... (RefTypes);
    }


    template<unsigned int Index>
    MemberReference<Class, type<Index>, refType<Index>>const getReference() const {
        return refs.MemRefs<Index, void>::ref;
    }

    template<unsigned int Index>
    type<Index> getValue(const Class* c) const {
        return getReference<Index>()->getValue(c);
    }

    template<unsigned int Index>
    type<Index> getValue(const Class& c) const {
        return getReference<Index>().getValue(c);
    }

    template<unsigned int Index>
    void setValue(Class* c, const type<Index>& val) const {
        getReference<Index>().setValue(c, val);
    }

    template<unsigned int Index>
    void setValue(Class& c, const type<Index>& val) const {
        getReference<Index>().setValue(c, val);
    }


    template<unsigned int Index>
    type<Index>& getAttr(Class* c) const {
        static_assert (IsDirectReference<memRefType<Index>>::value, "The current reference to the member does not provide direct approach.");
        return getReference<Index>()->getAttr(c);
    }

    template<unsigned int Index>
    type<Index>& getAttr(Class& c) const {
        static_assert (IsDirectReference<memRefType<Index>>::value, "The current reference to the member does not provide direct approach.");
        return getReference<Index>().getAttr(c);
    }


    template<unsigned int Index>
    const char* getName() const {
        return refs.MemRefs<Index, void>::name;
    }


    template<typename...Refs>
    Traits(Refs... refsAndNames) : refs(refsAndNames...){}
private:

    template<unsigned int Index = 0, typename Dummy = void>
    struct Apply {
        using ThisTrait = Traits<Class, RefTypes...>;

        template <class Functor>
        static auto apply (const ThisTrait& trait, Functor f,...)
        -> typename std::enable_if<std::is_assignable<
        std::function<
            void(unsigned int,
                 memRefType<Index>&,
                 const std::string&
                 )
        >, Functor>::value>::type
        {

            static_assert (std::is_assignable<
                    std::function<
                        void(unsigned int,
                             memRefType<Index>&,
                             const std::string&
                             )
                    >, Functor>::value, "");

            f(Index, trait.template getReference<Index>(), trait.template getName<Index>());
            Apply<Index + 1>::apply(trait, f);
        }



        template <class Functor>
        static auto apply (const ThisTrait& trait, Functor f)
        -> typename std::enable_if<std::is_assignable<
        std::function<
            void(
                 memRefType<Index>&,
                 const std::string&
                 )
        >, Functor>::value>::type
        {

            static_assert (std::is_assignable<
                    std::function<
                        void(
                             memRefType<Index>&,
                             const std::string&
                             )
                    >, Functor>::value, "");

            f(trait.template getReference<Index>(), trait.template getName<Index>());
            Apply<Index + 1>::apply(trait, f);
        }

        template <template <typename, typename>class Functor>
        static auto apply (const ThisTrait& trait)
        -> typename std::enable_if<std::is_class<Functor<Class, typename ThisTrait::template type<Index>>>::value>::type
        {

            static_assert (std::is_assignable<
                    std::function<
                        void(unsigned int,
                             memRefType<Index>&,
                             const std::string&
                             )
                    >, Functor<Class, typename ThisTrait::template type<Index>>>::value, "");


            Functor<Class, typename ThisTrait::template type<Index>>()(Index, trait.template getReference<Index>(), trait.template getName<Index>());
            Apply<Index + 1>::template apply<Functor>(trait);
        }

    };

    template<typename Dummy>
    struct Apply<size() - 1, Dummy> {
        using ThisTrait = Traits<Class, RefTypes...>;

        template <class Functor>
        static auto apply (const ThisTrait& trait, Functor f,...)
        -> typename std::enable_if<std::is_assignable<
        std::function<
            void(unsigned int,
                 memRefType<size() - 1>&,
                 const std::string&
                 )
        >, Functor>::value>::type
        {

            static_assert (std::is_assignable<
                    std::function<
                        void(unsigned int,
                             memRefType<size() - 1>&,
                             const std::string&
                             )
                    >, Functor>::value, "");

            f(trait.size() - 1, trait.template getReference<ThisTrait::size() - 1>(), trait.template getName<ThisTrait::size() - 1>());

        }



        template <class Functor>
        static auto apply (const ThisTrait& trait, Functor f)
        -> typename std::enable_if<std::is_assignable<
        std::function<
            void(
                 memRefType<size() - 1>&,
                 const std::string&
                 )
        >, Functor>::value>::type
        {

            static_assert (std::is_assignable<
                    std::function<
                        void(
                             memRefType<size() - 1>&,
                             const std::string&
                             )
                    >, Functor>::value, "");

            f(trait.template getReference<ThisTrait::size() - 1>(), trait.template getName<ThisTrait::size() - 1>());

        }

        template <template <typename, typename>class Functor>
        static auto apply (const ThisTrait& trait)
        -> typename std::enable_if<std::is_class<Functor<Class, typename ThisTrait::template type<ThisTrait::size() - 1>>>::value>::type
        {

            static_assert (std::is_assignable<
                    std::function<
                        void(unsigned int,
                             memRefType<size() - 1>&,
                             const std::string&
                             )
                    >, Functor<Class, typename ThisTrait::template type<ThisTrait::size() - 1>>>::value, "");


            Functor<Class, typename ThisTrait::template type<ThisTrait::size() - 1>>()(ThisTrait::size() - 1, trait.template getReference<ThisTrait::size() - 1>(), trait.template getName<ThisTrait::size() - 1>());
        }

    };
public:
    /**
     * @brief Traits apply function.
     * This function automatically
     * applies a lambda with specified
     * arguments: <BR> (unsigned int,
     * const auto& [as const MemberApproach<Class, typename>&]
     * const std::string&)
     */
    template<typename Functor>
        void apply(Functor f) const {
            Apply<>::apply(*this,f);
        }

    template<template <typename, typename>class Functor>
        void apply() const {
            Apply<>::template apply<Functor>(*this);
        }
};







template<typename Class>
class Traits<Class>{
    //static_assert (false, "The Traits template must be specialized for given type and must contain Traits references using variadic Traits.");
public:
    static constexpr std::false_type is_specialized{};
};




template<typename Class>
class DefaultIOTraits : public Traits<Class> {};

template<typename Class>
class DefaultArithmeticTraits : public Traits<Class> {};


#include "CustomTypeTraits.h"

namespace Impl {
template <unsigned int Index, unsigned int ...Indexes>
struct TraitedAttributeGetter{
    template<typename TraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value>::type>
    static auto& get(TraitT& arg){
        return TraitedAttributeGetter<Indexes...>::get(DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(arg));
    }

    template<typename TraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value>::type>
    static auto& get(TraitT* arg){
        return TraitedAttributeGetter<Indexes...>::get(DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(arg));
    }
};

template <unsigned int Index>
struct TraitedAttributeGetter<Index>{
    template<typename TraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value>::type>
    static auto& get(TraitT& arg){
        return DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(arg);
    }

    template<typename TraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value>::type>
    static auto& get(TraitT* arg){
        return DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(arg);
    }
};

} //Impl


template <unsigned int ...Indexes, typename ArythmeticTraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<ArythmeticTraitT>::value>::type>
auto& getTraitedAttribute(ArythmeticTraitT& arg){
    return Impl::TraitedAttributeGetter<Indexes...>::get(arg);
}


template <unsigned int ...Indexes, typename ArythmeticTraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<ArythmeticTraitT>::value>::type>
auto& getTraitedAttribute(ArythmeticTraitT* arg){
    return Impl::TraitedAttributeGetter<Indexes...>::get(arg);
}





#include "../Macros/MacroForEach.h"

#define IMPL_MEMREF_TYPE_CUSTOM(name, memberRef) decltype(memberRef)
#define IMPL_NAME_AND_REF(Class, name, member) name, &Class::member
#define IMPL_NAME_ATT(attribute) #attribute, attribute

#define IMPL_MAKE_CUSTOM_ATTRIBUTE_TRAIT(TraitName,Class,...) \
template<> \
class TraitName<Class>{ \
public: \
    static constexpr std::true_type is_specialized{}; \
    using traitsType = ::Traits<Class, FOR_EACH_2ARGS(IMPL_MEMREF_TYPE_CUSTOM, __VA_ARGS__)>; \
    static const traitsType getTraits() {return traitsType(__VA_ARGS__);} \
    static constexpr unsigned int size() {return traitsType::size();}\
};



#define MAKE_CUSTOM_ATTRIBUTE_TRAIT(Class,...) IMPL_MAKE_CUSTOM_ATTRIBUTE_TRAIT(Traits, Class, __VA_ARGS__) // defining specialization for Traits

#define MAKE_NAMED_ATTRIBUTE_TRAIT(Class, ...) MAKE_CUSTOM_ATTRIBUTE_TRAIT(Class, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TRAIT(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT(Class, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))


#define MAKE_CUSTOM_ATTRIBUTE_TRAIT_IO(Class,...) IMPL_MAKE_CUSTOM_ATTRIBUTE_TRAIT(DefaultIOTraits, Class,__VA_ARGS__) // defining specialization for DefaultIOTraits

#define MAKE_NAMED_ATTRIBUTE_TRAIT_IO(Class, ...) MAKE_CUSTOM_ATTRIBUTE_TRAIT_IO(Class, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TRAIT_IO(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT_IO(Class, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))


#define MAKE_CUSTOM_ATTRIBUTE_TRAIT_ARITHMETIC(Class,...) IMPL_MAKE_CUSTOM_ATTRIBUTE_TRAIT(DefaultArithmeticTraits, Class,__VA_ARGS__) // defining specialization for Traits

#define MAKE_NAMED_ATTRIBUTE_TRAIT_ARITHMETIC(Class, ...) MAKE_CUSTOM_ATTRIBUTE_TRAIT_ARITHMETIC(Class, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT_ARITHMETIC(Class, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))

#endif // TRAITS_H
