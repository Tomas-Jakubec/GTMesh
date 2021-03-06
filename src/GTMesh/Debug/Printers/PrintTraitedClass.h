#ifndef PRINTTRAITEDCLASS_H
#define PRINTTRAITEDCLASS_H
#include "../VariableExport.h"

struct PrintTraitedClass {

    template<typename T, typename TraitsIO, unsigned int Index = 0, bool = Index == TraitsIO::size() - 1>
    struct PrintClass{
        static void print(std::ostream& ost, const T &traitedClass, const TraitsIO& traitsIO){
            PrintClass<T, TraitsIO, Index, true>::print(ost, traitedClass, traitsIO);
            ost << ", ";
            PrintClass<T, TraitsIO, Index + 1>::print(ost, traitedClass, traitsIO);

        }
        template <typename ... _Traits, typename _T = T>
        static void print(std::ostream& ost, const _T &traitedClass, const TraitsIO& traitsIO, std::tuple<const _Traits& ...> t){
            PrintClass<T, TraitsIO, Index, true>::print(ost, traitedClass, traitsIO, t);
            ost << ", ";
            PrintClass<T, TraitsIO, Index + 1>::print(ost, traitedClass, traitsIO, t);

        }

        static void print(const T &traitedClass, const TraitsIO& traitsIO){
            PrintClass<T, TraitsIO, Index, true>::print(traitedClass, traitsIO);
            printf(", ");
            PrintClass<T, TraitsIO, Index + 1>::print(traitedClass, traitsIO);
        }

        template <typename ... _Traits, typename _T = T>
        static void print(const _T &traitedClass, const TraitsIO& traitsIO, std::tuple<const _Traits& ...> t){
            PrintClass<T, TraitsIO, Index, true>::print(traitedClass, traitsIO, t);
            printf(", ");
            PrintClass<T, TraitsIO, Index + 1>::print(traitedClass, traitsIO, t);
        }
    };

    template<typename T, typename TraitsIO, unsigned int Index>
    struct PrintClass<T, TraitsIO, Index, true>{

        static void print(std::ostream& ost, const T &traitedClass, const TraitsIO& traitsIO){
            ost << '"' << traitsIO.template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, traitsIO.template getValue<Index>(traitedClass));
        }

        template <typename ... _Traits, typename _T = T>
        static void print(std::ostream& ost, const _T &traitedClass, const TraitsIO& traitsIO, std::tuple<const _Traits& ...> t) {
            ost << '"' << traitsIO.template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, traitsIO.template getValue<Index>(traitedClass), t);
        }

        static void print(const T &traitedClass, const TraitsIO& traitsIO){
            printf("\"%s\" : ", traitsIO.template getName<Index>());
            VariableExport::exportVariable(traitsIO.template getValue<Index>(traitedClass));
        }

        template <typename ... _Traits, typename _T = T>
        static void print(const _T &traitedClass, const TraitsIO& traitsIO, std::tuple<const _Traits& ...> t) {
            printf("\"%s\" : ", traitsIO.template getName<Index>());
            VariableExport::exportVariable(traitsIO.template getValue<Index>(traitedClass), t);
        }
    };

    template<typename TraitedClass, std::enable_if_t<HasDefaultIOTraits<TraitedClass>::value, bool> = true>
    static void print(std::ostream& ost, const TraitedClass& var, ...)
    {
        ost << "{ ";
        PrintClass<TraitedClass, typename DefaultIOTraits<TraitedClass>::traitsType>::print(ost, var, DefaultIOTraits<TraitedClass>::getTraits());
        ost << " }";
    }

    template<typename Class, typename ...TraitsTypes>
    using IsValid = std::enable_if_t<SelectTraitsWithIODefault<Class, TraitsTypes...>::valid, bool>;


    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits>
    static void print(std::ostream& ost, const TraitsBinder<Class, PrimaryTraits, SecondaryTraits...> &traitedClass, ...)
    {
        VariableExport::exportVariable(ost, traitedClass.object, traitedClass.tupTraits);
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits, IsValid<Class, PrimaryTraits, SecondaryTraits...> = true>
    static void print(std::ostream& ost, const Class &traitedClass, const std::tuple<PrimaryTraits, SecondaryTraits...>& tupleTraits)
    {
        auto traits = SelectTraitsWithIODefault<Class, PrimaryTraits, SecondaryTraits...>::getTraitsInstance(tupleTraits);
        ost << "{ ";
        PrintClass<Class, decltype(traits)>::print(ost, traitedClass, traits, tupleTraits);
        ost << " }";
    }



    template<typename TraitedClass, std::enable_if_t<HasDefaultIOTraits<TraitedClass>::value, bool> = true>
    static void print(const TraitedClass& var, ...)
    {
        printf("{ ");
        PrintClass<TraitedClass, typename DefaultIOTraits<TraitedClass>::traitsType>::print(var, DefaultIOTraits<TraitedClass>::getTraits());
        printf(" }");
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits>
    static void print(const TraitsBinder<Class, PrimaryTraits, SecondaryTraits...> &traitedClass, ...)
    {
        VariableExport::exportVariable(traitedClass.object, traitedClass.tupTraits);
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits, IsValid<Class, PrimaryTraits, SecondaryTraits...> = true>
    static void print(const Class &traitedClass, const std::tuple<PrimaryTraits, SecondaryTraits...>& tupleTraits)
    {
        auto traits = SelectTraitsWithIODefault<Class, PrimaryTraits, SecondaryTraits...>::getTraitsInstance(tupleTraits);
        printf("{ ");
        PrintClass<Class, decltype(traits)>::print(traitedClass, traits, tupleTraits);
        printf(" }");
    }

};

#endif // PRINTTRAITEDCLASS_H
