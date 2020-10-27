#ifndef TRAITS_H
#define TRAITS_H
#include "MemberAccess/MemberAccess.h"
#include <string>
#include <memory>
#include <functional>

/**
 *
 */
template<typename Class, typename...RefTypes>
class Traits {
public:
    static constexpr bool isTraits = true;

    template<unsigned int Index>
    using refType = typename std::tuple_element<Index,std::tuple<RefTypes...>>::type;

    template<unsigned int Index>
    using memRefType = MemberAccess<Class, refType<Index>>;

    template <unsigned int Index>
    using type = typename MemberAccess<Class, refType<Index>>::typeValue;

private:

    template<unsigned int Index = 0, typename = void>
    struct MemRefs: public MemRefs<Index + 1> {

        const MemberAccess<Class, refType<Index>> ref;
        const char* name;

        template <typename ... REST>
        constexpr MemRefs(const char* n, refType<Index> r, REST... rest) : MemRefs<Index + 1> (rest...), ref(r), name(n){}
    };

    template<typename Dummy>
    struct MemRefs<sizeof...(RefTypes) - 1, Dummy>{
        const MemberAccess<Class, refType<sizeof...(RefTypes) - 1>> ref;
        const char* name;

        constexpr MemRefs(const char* n, refType<sizeof...(RefTypes) - 1> r) : ref(r), name(n){}
    };

    /**
     * @brief Container matching MemberAccess and a name.
     */
    const MemRefs<0, void> mReferences;

public:
    /**
     * The constructor of Traits initializes the names
     * and references in the refs container.
     * The first parameter is specified in order to distinguish the
     * copy and move constructors.
     * @param refsAndNames a parameter pack of names and references
     */
    template<typename...Refs>
    constexpr Traits(const char* n1, Refs... refsAndNames) : mReferences(n1, refsAndNames...){}

    constexpr Traits(const Traits<Class, RefTypes...>&) = default;
    constexpr Traits(Traits<Class, RefTypes...>&&) = default;

    static constexpr unsigned int size(){
        return sizeof... (RefTypes);
    }


    template<unsigned int Index>
    inline const MemberAccess<Class, refType<Index>>& getReference() const {
            return mReferences.MemRefs<Index, void>::ref;
    }

    template<unsigned int Index>
    inline type<Index> getValue(Class* c) const {
        return getReference<Index>().getValue(c);
    }

    template<unsigned int Index>
    inline type<Index> getValue(Class& c) const {
        return getReference<Index>().getValue(c);
    }

    template<unsigned int Index>
    inline type<Index> getValue(const Class* c) const {
        static_assert (HasConstGetAccess<memRefType<Index>>::value, "The current reference to the member does not provide constant approach.");
        return getReference<Index>().getValue(c);
    }

    template<unsigned int Index>
    inline type<Index> getValue(const Class& c) const {
        static_assert (HasConstGetAccess<memRefType<Index>>::value, "The current reference to the member does not provide constant approach.");
        return getReference<Index>().getValue(c);
    }

    template<unsigned int Index>
    inline void setValue(Class* c, const type<Index>& val) const {
        getReference<Index>().setValue(c, val);
    }

    template<unsigned int Index>
    inline void setValue(Class& c, const type<Index>& val) const {
        getReference<Index>().setValue(c, val);
    }


    template<unsigned int Index>
    inline type<Index>& getAttr(Class* c) const {
        static_assert (IsDirectAccess<memRefType<Index>>::value, "The current reference to the member does not provide direct approach.");
        return getReference<Index>().getAttr(c);
    }

    template<unsigned int Index>
    inline type<Index>& getAttr(Class& c) const {
        static_assert (IsDirectAccess<memRefType<Index>>::value, "The current reference to the member does not provide direct approach.");
        return getReference<Index>().getAttr(c);
    }


    template<unsigned int Index>
    inline const char* getName() const {
        return mReferences.MemRefs<Index, void>::name;
    }


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

namespace Impl {

template<typename,typename >
struct MakeTraitsFromTuple {};

template<typename TraitedClass, typename ... Refs>
struct MakeTraitsFromTuple<TraitedClass, std::tuple<Refs...>> {
    using type = Traits<TraitedClass, Refs...>;
};

template<typename T1, typename A>
struct TupleAppend {};

template<typename T1, typename ... A>
struct TupleAppend<T1, std::tuple<A...>> {
    using type = std::tuple<T1, A...>;
};

template<typename odd, typename even, typename ... Refs>
struct SelectEvenTypes {
    using type = typename TupleAppend<even, typename SelectEvenTypes<Refs...>::type>::type;
};

template<typename odd, typename even>
struct SelectEvenTypes<odd, even> {
    using type = std::tuple<even>;
};

}

template<typename TraitedClass, typename ... Args>
auto
makeTraits(const Args&... args) {
    return typename Impl::MakeTraitsFromTuple<TraitedClass, typename Impl::SelectEvenTypes<Args...>::type>::type(args...);
}


template<typename Class>
class DefaultTraits{};


template<typename Class>
class DefaultIOTraits : public DefaultTraits<Class> {};

template<typename Class>
class DefaultArithmeticTraits : public DefaultTraits<Class> {};


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
auto& get(ArythmeticTraitT& arg){
    return Impl::TraitedAttributeGetter<Indexes...>::get(arg);
}


template <unsigned int ...Indexes, typename ArythmeticTraitT, typename = typename std::enable_if<HasDefaultArithmeticTraits<ArythmeticTraitT>::value>::type>
auto& get(ArythmeticTraitT* arg){
    return Impl::TraitedAttributeGetter<Indexes...>::get(arg);
}





#include "../Macros/MacroForEach.h"

#define IMPL_MEMREF_TYPE_CUSTOM(name, memberRef) decltype(memberRef)
#define IMPL_NAME_AND_REF(Class, name, member) name, (&Class::member)
#define IMPL_NAME_ATT(attribute) #attribute, attribute

#define IMPL_NAME_AND_REF_TEMPLATE(Class, name, member) name, (&UNWRAP(Class)::member)

#define IMPL_MAKE_CUSTOM_TRAIT(TraitName,Class,...) \
template<> \
class TraitName<Class>{ \
public: \
    using traitsType = ::Traits<Class, FOR_EACH_2ARGS(IMPL_MEMREF_TYPE_CUSTOM, __VA_ARGS__)>; \
    static const traitsType getTraits() {return traitsType(__VA_ARGS__);} \
    static constexpr unsigned int size() {return traitsType::size();}\
}

#define IMPL_MAKE_CUSTOM_TEMPLATE_TRAIT(TraitName, TemplateParameters, Class,...) \
template<UNWRAP(TemplateParameters)> \
class TraitName<UNWRAP(Class)>{ \
public: \
    using traitsType = ::Traits<UNWRAP(Class), FOR_EACH_2ARGS(IMPL_MEMREF_TYPE_CUSTOM, __VA_ARGS__)>; \
    static const traitsType getTraits() {return traitsType(__VA_ARGS__);} \
    static constexpr unsigned int size() {return traitsType::size();}\
}

#define MAKE_CUSTOM_TRAIT(Class,...) IMPL_MAKE_CUSTOM_TRAIT(DefaultTraits, Class, __VA_ARGS__) // defining specialization for Traits

#define MAKE_NAMED_ATTRIBUTE_TRAIT(Class, ...) MAKE_CUSTOM_TRAIT(Class, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TRAIT(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT(Class, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))


#define MAKE_CUSTOM_TEMPLATE_TRAIT(Class, TemplateParameters, ...) IMPL_MAKE_CUSTOM_TEMPLATE_TRAIT(DefaultTraits, TemplateParameters, Class, __VA_ARGS__) // defining specialization for Traits

#define MAKE_NAMED_ATTRIBUTE_TEMPLATE_TRAIT(Class, TemplateParameters, ...) MAKE_CUSTOM_TEMPLATE_TRAIT(Class, TemplateParameters, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF_TEMPLATE, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TEMPLATE_TRAIT(Class, TemplateParameters, ...) MAKE_NAMED_ATTRIBUTE_TEMPLATE_TRAIT(Class, TemplateParameters, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))


#define MAKE_CUSTOM_TRAIT_IO(Class,...) IMPL_MAKE_CUSTOM_TRAIT(DefaultIOTraits, Class,__VA_ARGS__) // defining specialization for DefaultIOTraits

#define MAKE_NAMED_ATTRIBUTE_TRAIT_IO(Class, ...) MAKE_CUSTOM_TRAIT_IO(Class, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TRAIT_IO(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT_IO(Class, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))


#define MAKE_CUSTOM_TEMPLATE_TRAIT_IO(Class, TemplateParameters, ...) IMPL_MAKE_CUSTOM_TEMPLATE_TRAIT(DefaultIOTraits, TemplateParameters, Class, __VA_ARGS__) // defining specialization for Traits

#define MAKE_NAMED_ATTRIBUTE_TEMPLATE_TRAIT_IO(Class, TemplateParameters, ...) MAKE_CUSTOM_TEMPLATE_TRAIT_IO(Class, TemplateParameters, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF_TEMPLATE, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TEMPLATE_TRAIT_IO(Class, TemplateParameters, ...) MAKE_NAMED_ATTRIBUTE_TEMPLATE_TRAIT_IO(Class, TemplateParameters, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))



#define MAKE_CUSTOM_TRAIT_ARITHMETIC(Class,...) IMPL_MAKE_CUSTOM_TRAIT(DefaultArithmeticTraits, Class,__VA_ARGS__) // defining specialization for DefaultArithmeticTraits

#define MAKE_NAMED_ATTRIBUTE_TRAIT_ARITHMETIC(Class, ...) MAKE_CUSTOM_TRAIT_ARITHMETIC(Class, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(Class, ...) MAKE_NAMED_ATTRIBUTE_TRAIT_ARITHMETIC(Class, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))


#define MAKE_CUSTOM_TEMPLATE_TRAIT_ARITHMETIC(Class, TemplateParameters, ...) IMPL_MAKE_CUSTOM_TEMPLATE_TRAIT(DefaultArithmeticTraits, TemplateParameters, Class, __VA_ARGS__) // defining specialization for Traits

#define MAKE_NAMED_ATTRIBUTE_TEMPLATE_TRAIT_ARITHMETIC(Class, TemplateParameters, ...) MAKE_CUSTOM_TEMPLATE_TRAIT_ARITHMETIC(Class, TemplateParameters, FOR_EACH_3ARGS_1STAT(IMPL_NAME_AND_REF_TEMPLATE, Class, __VA_ARGS__))

#define MAKE_ATTRIBUTE_TEMPLATE_TRAIT_ARITHMETIC(Class, TemplateParameters, ...) MAKE_NAMED_ATTRIBUTE_TEMPLATE_TRAIT_ARITHMETIC(Class, TemplateParameters, FOR_EACH(IMPL_NAME_ATT, __VA_ARGS__))

#endif // TRAITS_H
