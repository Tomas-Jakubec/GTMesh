#ifndef SERIALIZETRAITEDCLASS_H
#define SERIALIZETRAITEDCLASS_H
#include <GTMesh/Utils/ConstexprFor.h>
#include <GTMesh/Traits/Traits.h>
#include <GTMesh/Traits/TraitsBind/TraitsBind.h>
#include "../BinarySerializer.h"


struct SerializeTraitedClass{
protected:
    // Helper functions serializing the instance of traited class
    struct TraitedClassBinarySerializer
    {
        template<size_t Index, typename TraitedType, typename TraitsType>
        static void exec(BinarySerializer::ByteContainer &container,
                         const TraitedType &var,
                         const TraitsType &traits)
        {
            BinarySerializer::addNext(container, traits.template getValue<Index>(var));
        }

        template<size_t Index, typename TraitedType, typename TraitsType, typename ...TupleTraits>
        static void exec(BinarySerializer::ByteContainer &container,
                         const TraitedType &var,
                         const TraitsType &traits,
                         const std::tuple<TupleTraits...>& tupleTraits)
        {
            BinarySerializer::addNext(container, traits.template getValue<Index>(var));
        }
    };

    struct TraitedClassBinaryDeserializer
    {
        template<size_t Index, typename TraitedType, typename TraitsType>
        static void exec(BinarySerializer::ByteContainerIterator &it,
                         TraitedType &var,
                         const TraitsType &traits)
        {
            BinarySerializer::readNext(it, traits.template getAttr<Index>(var));
        }

        template<size_t Index, typename TraitedType, typename TraitsType,  typename ...TupleTraits>
        static void exec(BinarySerializer::ByteContainerIterator &it,
                         TraitedType &var,
                         const TraitsType &traits,
                         const std::tuple<TupleTraits...>& tupleTraits)
        {
            BinarySerializer::readNext(it, traits.template getAttr<Index>(var));
        }
    };

public:
    template<typename T,
             typename...,
             std::enable_if_t<HasDefaultArithmeticTraits<T>::value, bool> = true>
    static void serialize(BinarySerializer::ByteContainer &container, const T &var)
    {
        constexprFor<TraitedClassBinarySerializer, DefaultArithmeticTraits<T>::size()>(
            container, var, DefaultArithmeticTraits<T>::getTraits());
    }

    template<typename T,
             typename ...TupleTraits,
             std::enable_if_t<SelectTraitsWithArithmeticDefault<T, TupleTraits...>::valid, bool> = true>
    static void serialize(BinarySerializer::ByteContainer &container, const T &var,
                          const std::tuple<TupleTraits...>& tupleTraits)
    {
        constexprFor<TraitedClassBinarySerializer,
                     SelectTraitsWithArithmeticDefault<T, TupleTraits...>::TypeTraits::size()>(
            container,
            var,
            SelectTraitsWithArithmeticDefault<T, TupleTraits...>::getTraitsInstance(tupleTraits),
            tupleTraits);
    }

    template<typename T,
             typename PrimaryTraits,
             typename ...SecondaryTraits>
    static void serialize(BinarySerializer::ByteContainer &container, const TraitsBinder<T, PrimaryTraits, SecondaryTraits...> &traitedClass)
    {
        BinarySerializer::addNext(container, traitedClass.object, traitedClass.tupTraits);
    }

    template<typename T,
             typename...,
             std::enable_if_t<HasDefaultArithmeticTraits<T>::value, bool> = true>
    static void deserialize(BinarySerializer::ByteContainerIterator &it, T &var)
    {
        constexprFor<TraitedClassBinaryDeserializer, DefaultArithmeticTraits<T>::size()>(
            it, var, DefaultArithmeticTraits<T>::getTraits());
    }

    template<typename T,
             typename ...TupleTraits,
             std::enable_if_t<SelectTraitsWithArithmeticDefault<T, TupleTraits...>::valid, bool> = true>
    static void deserialize(BinarySerializer::ByteContainerIterator &it, T &var,
                            const std::tuple<TupleTraits...>& tupleTraits)
    {
        constexprFor<TraitedClassBinaryDeserializer,
                     SelectTraitsWithArithmeticDefault<T, TupleTraits...>::TypeTraits::size()>(
            it,
            var,
            SelectTraitsWithArithmeticDefault<T, TupleTraits...>::getTraitsInstance(tupleTraits),
            tupleTraits);
    }

    template<typename T,
             typename PrimaryTraits,
             typename ...SecondaryTraits>
    static void deserialize(BinarySerializer::ByteContainerIterator &it, const TraitsBinder<T, PrimaryTraits, SecondaryTraits...> &traitedClass)
    {
        BinarySerializer::readNext(it, traitedClass.object, traitedClass.tupTraits);
    }
};
#endif // SERIALIZETRAITEDCLASS_H
