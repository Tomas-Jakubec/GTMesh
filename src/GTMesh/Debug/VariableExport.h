#ifndef VARIABLEEXPORT_H
#define VARIABLEEXPORT_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../Traits/Traits.h"
#include "../Traits/CustomTypeTraits.h"


enum VARIABLE_EXPORT_METHOD {
    ostream,
    stdio
};

template <VARIABLE_EXPORT_METHOD target = VARIABLE_EXPORT_METHOD::ostream>
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
            !std::is_same<T, char>::value &&
            !HasDefaultIOTraits<T>::value
         >::type
    {
        ost << b;
    }



    static void exportVariable(std::ostream& ost, const bool& b)
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
    static auto exportVariable(std::ostream& ost, const T &list)
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
    static auto exportVariable(std::ostream& ost, const T &list)
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


    template<typename T, unsigned int Index = 0, bool = Index == DefaultIOTraits<T>::size() - 1>
    struct PrintClass{
        static void print(std::ostream& ost, const T &traitedClass){
            PrintClass<T, Index, true>::print(ost, traitedClass);
            ost << ", ";
            PrintClass<T, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T, unsigned int Index>
    struct PrintClass<T, Index, true>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << DefaultIOTraits<T>::getTraits().template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, DefaultIOTraits<T>::getTraits().template getValue<Index>(traitedClass));
        }
    };


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &traitedClass)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        ost << "{ ";
        PrintClass<T>::print(ost, traitedClass);
        ost << " }";
    }


};


template <>
struct VariableExport<VARIABLE_EXPORT_METHOD::stdio> {


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
    static auto exportVariable(const std::pair<T1,T2>& b) -> void
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


    template<typename T, unsigned int Index = 0, bool = Index == DefaultIOTraits<T>::size() - 1>
    struct PrintClass{
        static void print(const T &traitedClass){
            PrintClass<T, Index, true>::print(traitedClass);
            printf(", ");
            PrintClass<T, Index + 1>::print(traitedClass);

        }
    };

    template<typename T, unsigned int Index>
    struct PrintClass<T, Index, true>{
        static void print(const T &traitedClass){
            printf("\"%s\" : ", DefaultIOTraits<T>::getTraits().template getName<Index>());
            VariableExport::exportVariable(DefaultIOTraits<T>::getTraits().template getValue<Index>(traitedClass));
        }
    };



    template<typename T>
    static auto exportVariable(const T &traitedClass)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        printf("{ ");
        PrintClass<T>::print(traitedClass);
        printf(" }");
    }


};

#endif // VARIABLEEXPORT_H