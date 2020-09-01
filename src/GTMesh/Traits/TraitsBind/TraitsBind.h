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

#endif // TRAITSBIND_H
