#ifndef VARIABLEEXPORT_H
#define VARIABLEEXPORT_H
#include <iostream>
#include <fstream>
#include <string>
#include "../Traits/Traits.h"
#include "../Traits/CustomTypeTraits.h"

struct VariableExport {


    static void exportVariable(std::ostream& ost, ...)
    {
        ost << "\"variable is not exportable\"" << std::endl;
    }


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T& b)
      -> typename std::enable_if<
             IsExportable<T>::value &&
            !std::is_same<T, bool>::value &&
            !std::is_same<T, std::string>::value &&
            !std::is_same<T, const char*>::value &&
            !std::is_same<T, char*>::value &&
            !std::is_same<T, const char>::value
         >::type
    {
        ost << b;
    }



    static void exportVariable(std::ostream& ost, const bool& b)
    {
        ost << (b == true ? "true" : "false");
    }

    static void exportVariable(std::ostream& ost, const std::string& str)
    {
        ost << '"' << str << '"';
    }


    static void exportVariable(std::ostream& ost, const char* str)
    {
        ost << '"' << str << '"';
    }


    static void exportVariable(std::ostream& ost, const char str)
    {
        ost << '"' << str << '"';
    }


    template<typename T1, typename T2>
    static auto exportVariable(std::ostream& ost, const std::pair<T1,T2>& b) -> void
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
             !HasDefaultTraits<T>::value
         >::type
    {
        auto it = list.begin();
        ost << "[ ";
        while (it != list.end()){
            exportVariable(ost, *it);
            if (++it != list.end()){
                ost << ", ";
            }
        }
        ost << " ]";
    }



    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultTraits<T>::value
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
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsTNLIndexable<T>::value &&
             !IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultTraits<T>::value
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
    static void exportVariable(std::ostream& ost, const std::initializer_list<T> &list)
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


    template<typename T,unsigned int Index = 0, typename Void = void>
    struct PrintClass{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T>::ttype::template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, Traits<T>::ttype::template getReference<Index>()->getValue(traitedClass));
            ost << ", ";
            PrintClass<T, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T,unsigned int Index, typename... Types>
    struct PrintClass <Traits<T, Types...>, Index, typename std::enable_if<Index < Traits<T, Types...>::size() - 1>::type>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T, Types...>::template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, Traits<T, Types...>::template getReference<Index>()->getValue(traitedClass));
            ost << ", ";
            PrintClass<Traits<T, Types...>, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct PrintClass <Traits<T, Types...>, Index, typename std::enable_if<Index == Traits<T, Types...>::size() - 1>::type>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T, Types...>::template getName<Traits<T, Types...>::size() - 1>() << "\" : ";
            VariableExport::exportVariable(ost, Traits<T, Types...>::template getReference<Traits<T, Types...>::size() - 1>()->getValue(traitedClass));
        }
    };

    template<typename T, unsigned int Index>
    struct PrintClass<T, Index, typename std::enable_if<Index == Traits<T>::ttype::size() - 1>::type>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T>::ttype::template getName<Traits<T>::ttype::size() - 1>() << "\" : ";
            VariableExport::exportVariable(ost, Traits<T>::ttype::template getReference<Traits<T>::ttype::size() - 1>()->getValue(traitedClass));
        }
    };


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &traitedClass)
      -> typename std::enable_if<
             HasDefaultTraits<T>::value
         >::type
    {
        ost << "{ ";
        PrintClass<typename Traits<T>::ttype>::print(ost, traitedClass);
        ost << " }";
    }


};

#endif // VARIABLEEXPORT_H
