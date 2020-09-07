#ifndef PRINTTRAITEDCLASS_H
#define PRINTTRAITEDCLASS_H
#include "../VariableExport.h"

template <typename Class, typename ClassTraits>
struct IsTraitsOf : public std::false_type {};

template <typename Class, typename ClassTraits>
struct IsTraitsOf<Class, const ClassTraits&> : public IsTraitsOf<Class, std::decay_t<ClassTraits>>{};


template <typename Class, typename ... TraitsArgs>
struct IsTraitsOf<Class, Traits<Class, TraitsArgs...>> : public std::true_type {};


/**
 * This class selects the first class traits with respect
 * to Class.
 */
template <typename Class, size_t Index, typename TraitsType, typename... TraitsTypes>
struct SelectTraits : public
        std::conditional_t<IsTraitsOf<Class, TraitsType>::value, SelectTraits<Class, Index, TraitsType>, SelectTraits<Class, Index + 1, TraitsTypes...>>{};

template <typename Class, size_t Index, typename TraitsType>
struct SelectTraits<Class, Index, TraitsType> {
    static constexpr bool valid = IsTraitsOf<Class, TraitsType>::value || HasDefaultIOTraits<Class>::value;
private:
    template <bool b, bool valid, typename = void>
    struct _conditional{
        using type = TraitsType;
    };

    template <typename Dummy>
    struct _conditional<false, true, Dummy>{
        using type = typename DefaultIOTraits<Class>::traitsType;
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
        return DefaultIOTraits<Class>::getTraits();
    }
};

struct PrintTraitedClass {
    static int print(...){return 0;}

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
            VariableExport<>::exportVariable(ost, traitsIO.template getValue<Index>(traitedClass));
        }

        template <typename ... _Traits, typename _T = T>
        static void print(std::ostream& ost, const _T &traitedClass, const TraitsIO& traitsIO, std::tuple<const _Traits& ...> t) {
            ost << '"' << traitsIO.template getName<Index>() << "\" : ";
            VariableExport<>::exportVariable(ost, traitsIO.template getValue<Index>(traitedClass), t);
        }

        static void print(const T &traitedClass, const TraitsIO& traitsIO){
            printf("\"%s\" : ", traitsIO.template getName<Index>());
            VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(traitsIO.template getValue<Index>(traitedClass));
        }

        template <typename ... _Traits, typename _T = T>
        static void print(const _T &traitedClass, const TraitsIO& traitsIO, std::tuple<const _Traits& ...> t) {
            printf("\"%s\" : ", traitsIO.template getName<Index>());
            VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(traitsIO.template getValue<Index>(traitedClass), t);
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
    using IsValid = std::enable_if_t<SelectTraits<Class, 0, TraitsTypes...>::valid, bool>;


    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits>
    static void print(std::ostream& ost, const TraitsBinder<Class, PrimaryTraits, SecondaryTraits...> &traitedClass, ...)
    {
        VariableExport<>::exportVariable(ost, traitedClass.object, traitedClass.tupTraits);
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits, IsValid<Class, PrimaryTraits, SecondaryTraits...> = true>
    static void print(std::ostream& ost, const Class &traitedClass, const std::tuple<PrimaryTraits, SecondaryTraits...>& tupleTraits)
    {
        auto traits = SelectTraits<Class, 0, PrimaryTraits, SecondaryTraits...>::getTraitsInstance(tupleTraits);
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
        VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(traitedClass.object, traitedClass.tupTraits);
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits, IsValid<Class, PrimaryTraits, SecondaryTraits...> = true>
    static void print(const Class &traitedClass, const std::tuple<PrimaryTraits, SecondaryTraits...>& tupleTraits)
    {
        auto traits = SelectTraits<Class, 0, PrimaryTraits, SecondaryTraits...>::getTraitsInstance(tupleTraits);
        printf("{ ");
        PrintClass<Class, decltype(traits)>::print(traitedClass, traits, tupleTraits);
        printf(" }");
    }
};

#endif // PRINTTRAITEDCLASS_H
