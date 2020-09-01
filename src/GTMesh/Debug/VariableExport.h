#ifndef VARIABLEEXPORT_H
#define VARIABLEEXPORT_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../Traits/Traits.h"
#include "../Traits/TraitsBind/TraitsBind.h"
#include "../Traits/CustomTypeTraits.h"


template <typename Class, typename ClassTraits>
struct isTraitsOf : public std::false_type {};

template <typename Class, typename ClassTraits>
struct isTraitsOf<Class, const ClassTraits&> : public isTraitsOf<Class, std::decay_t<ClassTraits>>{};


template <typename Class, typename ... TraitsArgs>
struct isTraitsOf<Class, Traits<Class, TraitsArgs...>> : public std::true_type {};



template <typename Class,typename TraitsType, typename... TraitsTypes>
struct SelectTraits : public
        std::conditional_t<isTraitsOf<Class, TraitsType>::value, SelectTraits<Class, TraitsType>, SelectTraits<Class, TraitsTypes...>>{};


template <typename Class,typename TraitsType>
struct SelectTraits<Class, TraitsType> {
    using TypeTraits = std::conditional_t< isTraitsOf<Class, TraitsType>::value, TraitsType, typename DefaultIOTraits<Class>::traitsType>;

    template<typename ... Args, typename TT = TraitsType, typename std::enable_if<isTraitsOf<Class, TT>::value, bool>::type = true>
    static const auto& getTraitsInstance(const std::tuple<Args...>& t) {
        return std::get<const TypeTraits&>(t);
    }


    template<typename ... Args, typename TT = TraitsType, typename std::enable_if<!isTraitsOf<Class, TT>::value, bool>::type = true>
    static auto getTraitsInstance(const std::tuple<Args...>&) {
        return DefaultIOTraits<Class>::getTraits();
    }
};



enum VARIABLE_EXPORT_METHOD {
    VARIABLE_EXPORT_METHOD_OSTREAM,
    VARIABLE_EXPORT_METHOD_STDIO
};

template <VARIABLE_EXPORT_METHOD target = VARIABLE_EXPORT_METHOD::VARIABLE_EXPORT_METHOD_OSTREAM>
struct VariableExport {


    static void exportVariable(std::ostream& ost, ...)
    {
        ost << "\"variable is not exportable\"" << std::endl;
    }


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T& b, ...)
      -> typename std::enable_if<
             IsExportable<T>::value &&
            !std::is_same<T, bool>::value &&
            !std::is_same<T, std::string>::value &&
            !std::is_same<T, const char*>::value &&
            !std::is_same<T, char*>::value &&
            !std::is_same<T, char>::value &&
            !HasDefaultIOTraits<T>::value
         >::type
    {
        ost << b;
    }



    static void exportVariable(std::ostream& ost, const bool& b, ...)
    {
        ost << (b == true ? "true" : "false");
    }

    template<typename T>
    static auto exportVariable(std::ostream& ost, const T& str)
    -> typename std::enable_if<
            std::is_same<T, std::string>::value ||
            std::is_same<T, const char*>::value ||
            std::is_same<T, char*>::value ||
            std::is_same<T, char>::value
        >::type
    {
        ost << '"' << str << '"';
    }


    template<typename T1, typename T2>
    static auto exportVariable(std::ostream& ost, const std::pair<T1,T2>& b, ...)
    -> typename std::enable_if<
            !IsTraits<T2>::value
       >::type
    {
        ost << "{ ";
        exportVariable(ost, b.first);
        ost << ": ";
        exportVariable(ost, b.second);
        ost << "}";
    }

    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultIOTraits<T>::value
         >::type
    {
        auto it = list.cbegin();
        ost << "[ ";
        while (it != list.cend()){
            exportVariable(ost, *it);
            if (++it != list.cend()){
                ost << ", ";
            }
        }
        ost << " ]";
    }



    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list, ...)
      -> typename std::enable_if<
              IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultIOTraits<T>::value
         >::type
    {
        ost << "[ ";
        for (decltype (list.size())i = 0; i < list.size(); i++){
            exportVariable(ost, list[i]);
            if (i < list.size() - 1){
                ost << ", ";
            }
        }
        ost << " ]";
    }


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list, ...)
      -> typename std::enable_if<
              IsTNLIndexable<T>::value &&
             !IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultIOTraits<T>::value
         >::type
    {
        ost << "[ ";
        for (decltype (list.size())i = 0; i < list.size(); i++){
            exportVariable(ost, list[i]);
            if (i <  list.getSize() - 1){
                ost << ", ";
            }
        }
        ost << " ]";
    }



    template<typename T>
    static void exportVariable(std::ostream& ost, const std::initializer_list<T> &list, ...)
    {
        static auto it = list.begin();
        ost << "[ ";
        while (it != list.end()){
            exportVariable(ost, *it);
            if (++it != list.end()){
                ost << ", ";
            }
        }
        ost << " ]";
    }


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
    };



    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &traitedClass, ...)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        ost << "{ ";
        PrintClass<T, typename DefaultIOTraits<T>::traitsType>::print(ost, traitedClass, DefaultIOTraits<T>::getTraits());
        ost << " }";
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits>
    static void exportVariable(std::ostream& ost, const TraitsBinder<Class, PrimaryTraits, SecondaryTraits...> &traitedClass)
    {
        auto traits = SelectTraits<Class, PrimaryTraits, SecondaryTraits...>::getTraitsInstance(traitedClass.tupTraits);
        ost << "{ ";
        PrintClass<Class, decltype(traits)>::print(ost, traitedClass.object, traits, traitedClass.tupTraits);
        ost << " }";
    }

    template< typename Class, typename PrimaryTraits,  typename ... SecondaryTraits, typename std::enable_if<std::is_class<Class>::value, bool>::type = true>
    static void exportVariable(std::ostream& ost, const Class &traitedClass, std::tuple<PrimaryTraits, SecondaryTraits...> tupleTraits)
    {
        auto traits = SelectTraits<Class, PrimaryTraits, SecondaryTraits...>::getTraitsInstance(tupleTraits);
        ost << "{ ";
        PrintClass<Class, decltype(traits)>::print(ost, traitedClass, traits, tupleTraits);
        ost << " }";
    }

    template<typename T, typename... Refs>
    static auto exportVariable(std::ostream& ost, const std::pair<T, Traits<std::decay_t<T>, Refs...>> &traitedClassWithTraits, ...)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        ost << "{ ";
        PrintClass<T, Traits<std::decay_t<T>, Refs...>>::print(ost, traitedClassWithTraits.first, traitedClassWithTraits.second);
        ost << " }";
    }

};


template <>
struct VariableExport<VARIABLE_EXPORT_METHOD::VARIABLE_EXPORT_METHOD_STDIO> {


    static void exportVariable(...)
    {
        printf("\"variable is not exportable\"");
    }

    static void exportVariable(const double& d)
    {
        printf("%g", d);
    }

    static void exportVariable(const long double& d)
    {
        printf("%Lg", d);
    }

    static void exportVariable(const int& d)
    {
        printf("%d", d);
    }


    static void exportVariable(const short& d)
    {
        printf("%hd", d);
    }

    static void exportVariable(const unsigned int& d)
    {
        printf("%ud", d);
    }


    static void exportVariable(const long long& d)
    {
        printf("%lld", d);
    }

    static void exportVariable(const unsigned long long& d)
    {
        printf("%llu", d);
    }


    static void exportVariable(const bool& b)
    {
        printf("%s", (b == true ? "true" : "false"));
    }

    static void exportVariable(const std::string& str)
    {
        printf("\"%s\"", str.c_str());
    }


    static void exportVariable(const char* str)
    {
        printf("\"%s\"", str);
    }


    static void exportVariable(const char str)
    {
        printf("\"%c\"", str);
    }


    template<typename T1, typename T2>
    static auto exportVariable(const std::pair<T1,T2>& b)
    -> typename std::enable_if<
            !IsTraits<T2>::value
       >::type
    {
        printf("{ ");
        exportVariable(b.first);
        printf(": ");
        exportVariable(b.second);
        printf("}");
    }

    template<typename T>
    static auto exportVariable(const T &list)
      -> typename std::enable_if<
              IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultIOTraits<T>::value
         >::type
    {
        auto it = list.begin();
        printf("[ ");
        while (it != list.end()){
            exportVariable(*it);
            if (++it != list.end()){
                printf(", ");
            }
        }
        printf(" ]");
    }



    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultIOTraits<T>::value
         >::type
    {
        printf("[ ");
        for (decltype (list.size())i = 0; i < list.size(); i++){
            exportVariable(ost, list[i]);
            if (i < list.size() - 1){
                printf(", ");
            }
        }
        printf(" ]");
    }


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsTNLIndexable<T>::value &&
             !IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultIOTraits<T>::value
         >::type
    {
        printf("[ ");
        for (decltype (list.size())i = 0; i < list.size(); i++){
            exportVariable(ost, list[i]);
            if (i <  list.getSize() - 1){
                printf(", ");
            }
        }
        printf(" ]");
    }



    template<typename T>
    static void exportVariable(std::ostream& ost, const std::initializer_list<T> &list)
    {
        static auto it = list.begin();
        printf("[ ");
        while (it != list.end()){
            exportVariable(ost, *it);
            if (++it != list.end()){
                printf(", ");
            }
        }
        printf(" ]");
    }


    template<typename T, typename TraitsIO, unsigned int Index = 0, bool = Index == TraitsIO::size() - 1>
    struct PrintClass{
        static void print(const T &traitedClass, const TraitsIO& traitsIO){
            PrintClass<T, TraitsIO, Index, true>::print(traitedClass, traitsIO);
            printf(", ");
            PrintClass<T, TraitsIO, Index + 1>::print(traitedClass, traitsIO);

        }
    };

    template<typename T, typename TraitsIO, unsigned int Index>
    struct PrintClass<T, TraitsIO, Index, true>{
        static void print(const T &traitedClass, const TraitsIO& traitsIO){
            printf("\"%s\" : ", traitsIO.template getName<Index>());
            VariableExport::exportVariable(traitsIO.template getValue<Index>(traitedClass));
        }
    };


    template<typename T>
    static auto exportVariable(const T &traitedClass)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        printf("{ ");
        PrintClass<T, typename DefaultIOTraits<T>::traitsType>::print(traitedClass, DefaultIOTraits<T>::getTraits());
        printf(" }");
    }

    template<typename T, typename... Refs>
    static auto exportVariable(const std::pair<T, Traits<std::decay_t<T>, Refs...>>& traitedClassWithTraits)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        printf("{ ");
        PrintClass<T, Traits<std::decay_t<T>, Refs...>>::print(traitedClassWithTraits.first, traitedClassWithTraits.second);
        printf(" }");
    }


};

#endif // VARIABLEEXPORT_H
