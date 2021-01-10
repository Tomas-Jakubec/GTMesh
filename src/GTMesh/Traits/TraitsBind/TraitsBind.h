#ifndef TRAITSBIND_H
#define TRAITSBIND_H
#include "../Traits.h"

/**
 * @brief Non-owning container binding an object to its Traits.
 * It is possible to setup more than one Traits in case of
 * deeper structure of the object, i.e., it has members with
 * another Traits assigned.
 */
template <typename Class, typename PrimaryTraits,  typename ... SecondaryTraits>
struct TraitsBinder {
    // The object the traits are bound to
    Class& object;
    // Custom Traits used for description of traited class object
    const std::tuple<const PrimaryTraits&, const SecondaryTraits& ...> tupTraits;


    // Constructors
    TraitsBinder(Class& obj, const PrimaryTraits& primaryTraits, const SecondaryTraits&... secondaryTraits)
        : object(obj), tupTraits(primaryTraits, secondaryTraits...)
    {}

    TraitsBinder(const TraitsBinder<Class, PrimaryTraits, SecondaryTraits...>&) = default;
    TraitsBinder(TraitsBinder<Class, PrimaryTraits, SecondaryTraits...>&&) = default;
};


/**
 * @brief This function creates an instance of TraitsBinder
 * which gathers an instance of a traited class and traits
 * to be used to annotate the data members.
 */
template <typename Class, typename PrimaryTraits,  typename ... SecondaryTraits>
TraitsBinder<Class, PrimaryTraits, SecondaryTraits ...> bindTraits(Class& obj, const PrimaryTraits& p, const SecondaryTraits& ... s) {
    return TraitsBinder<Class, PrimaryTraits, SecondaryTraits ...>(obj, p, s...);
}


// Selector of traits from binded traits
template <typename Class, typename ClassTraits>
struct IsTraitsOf : public std::false_type {};

template <typename Class, typename ClassTraits>
struct IsTraitsOf<Class, const ClassTraits&> : public IsTraitsOf<Class, std::decay_t<ClassTraits>>{};


template <typename Class, typename ... TraitsArgs>
struct IsTraitsOf<Class, Traits<Class, TraitsArgs...>> : public std::true_type {};


/**
 * This class selects the first class traits according
 * to Class. If none type of traits is given, it returns
 * DefaultIOTraits<Class> if possible.
 */
template <typename Class, size_t Index, template<typename> class DefaultTraitsType, typename... TraitsTypes>
struct SelectTraits{};

template<typename Class,
         size_t Index,
         template<typename>
         class DefaultTraitsType,
         typename TraitsType,
         typename... TraitsTypes>
struct SelectTraits<Class, Index, DefaultTraitsType, TraitsType, TraitsTypes...>
    : public std::conditional_t<IsTraitsOf<Class, TraitsType>::value,
                                SelectTraits<Class, Index, DefaultTraitsType, TraitsType>,
                                SelectTraits<Class, Index + 1, DefaultTraitsType, TraitsTypes...>>
{};

template <typename Class, size_t Index, template<typename> class DefaultTraitsType, typename TraitsType>
struct SelectTraits<Class, Index, DefaultTraitsType, TraitsType> {
    static constexpr bool valid = IsTraitsOf<Class, TraitsType>::value || HasDefaultTraitsOfType<Class, DefaultTraitsType>::value;
private:
    template <bool b, bool valid, typename = void>
    struct _conditional{
        using type = typename std::decay<TraitsType>::type;
    };

    template <typename Dummy>
    struct _conditional<false, true, Dummy>{
        using type = typename DefaultTraitsType<Class>::traitsType;
    };


    template <typename Dummy>
    struct _conditional<false, false, Dummy>{
        using type = void;
    };
public:
    using TypeTraits = typename _conditional< IsTraitsOf<Class, TraitsType>::value, valid >::type;

    template<typename ... Args, typename TT = TraitsType, typename std::enable_if<IsTraitsOf<Class, TT>::value, bool>::type = true>
    static const auto& getTraitsInstance(const std::tuple<Args...>& t) {
        return std::get<Index>(t);
    }


    template<typename ... Args, typename TT = TraitsType, typename std::enable_if<!IsTraitsOf<Class, TT>::value, bool>::type = true>
    static auto getTraitsInstance(const std::tuple<Args...>&) {
        return DefaultTraitsType<Class>::getTraits();
    }
};

template <typename Class, size_t Index, template<typename> class DefaultTraitsType>
struct SelectTraits<Class, Index, DefaultTraitsType> {
    static constexpr bool valid = HasDefaultTraitsOfType<Class, DefaultTraitsType>::value;
private:
    template <bool valid, typename = void>
    struct _conditional{
        using type = typename DefaultTraitsType<Class>::traitsType;
    };


    template <typename Dummy>
    struct _conditional<false, Dummy>{
        using type = void;
    };
public:
    using TypeTraits = typename _conditional< valid >::type;


    template<typename ... Args>
    static auto getTraitsInstance(const std::tuple<Args...>&) {
        return DefaultTraitsType<Class>::getTraits();
    }
};

template <typename Class, size_t Index, typename... TraitsTypes>
using SelectTraitsWithIODefault = SelectTraits<Class, Index, DefaultIOTraits, TraitsTypes...>;

template <typename Class, size_t Index, typename... TraitsTypes>
using SelectTraitsWithArithmeticDefault = SelectTraits<Class, Index, DefaultArithmeticTraits, TraitsTypes...>;


#endif // TRAITSBIND_H
