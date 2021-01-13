#ifndef SERIALIZETRAITEDCLASS_H
#define SERIALIZETRAITEDCLASS_H
#include <GTMesh/Utils/ConstexprFor.h>
#include <GTMesh/Traits/Traits.h>
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
             typename...,
             std::enable_if_t<HasDefaultArithmeticTraits<T>::value, bool> = true>
    static void deserialize(BinarySerializer::ByteContainerIterator &it, T &var)
    {
        constexprFor<TraitedClassBinaryDeserializer, DefaultArithmeticTraits<T>::size()>(
            it, var, DefaultArithmeticTraits<T>::getTraits());
    }
};
#endif // SERIALIZETRAITEDCLASS_H
