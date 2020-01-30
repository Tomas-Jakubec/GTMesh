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
            !std::is_same<T, const char>::value &&
            !HasDefaultIOTraits<T>::value
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
            ost << '"' << Traits<T>::tr.template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, Traits<T>::tr.template getValue<Index>(traitedClass));
            ost << ", ";
            PrintClass<T, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T,unsigned int Index, typename... Types>
    struct PrintClass <Traits<T, Types...>, Index, typename std::enable_if<Index < Traits<T, Types...>::size() - 1>::type>{
        static void print(std::ostream& ost, const T &traitedClass, const Traits<T, Types...>& trait = DefaultIOTraits<T>::getTraits()){
            ost << '"' << trait.template getName<Index>() << "\" : ";
            VariableExport::exportVariable(ost, trait.template getValue<Index>(traitedClass));
            ost << ", ";
            PrintClass<Traits<T, Types...>, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct PrintClass <Traits<T, Types...>, Index, typename std::enable_if<Index == Traits<T, Types...>::size() - 1>::type>{
        static void print(std::ostream& ost, const T &traitedClass, const Traits<T, Types...>& trait = DefaultIOTraits<T>::getTraits()){
            ost << '"' << trait.template getName<Traits<T, Types...>::size() - 1>() << "\" : ";
            VariableExport::exportVariable(ost, trait.template getValue<Traits<T, Types...>::size() - 1>(traitedClass));
        }
    };

    template<typename T, unsigned int Index>
    struct PrintClass<T, Index, typename std::enable_if<Index == Traits<T>::traitsType::size() - 1>::type>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T>::tr.template getName<Traits<T>::traitsType::size() - 1>() << "\" : ";
            VariableExport::exportVariable(ost, Traits<T>::tr.template getValue<Traits<T>::tr.size() - 1>(traitedClass));
        }
    };


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &traitedClass)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value
         >::type
    {
        ost << "{ ";
        PrintClass<typename DefaultIOTraits<T>::traitsType>::print(ost, traitedClass);
        ost << " }";
    }


};


template <>
struct VariableExport<VARIABLE_EXPORT_METHOD::stdio> {


    static std::string exportVariable(...)
    {
        return "\"variable is not exportable\"";
    }


    template<typename T>
    static auto exportVariable(const T& b)
      -> typename std::enable_if<
             IsExportable<T>::value &&
            !std::is_same<T, bool>::value &&
            !std::is_same<T, std::string>::value &&
            !std::is_same<T, const char*>::value &&
            !std::is_same<T, char*>::value &&
            !std::is_same<T, const char>::value,
            std::string
         >::type
    {
        std::stringstream ss;
        ss << b;
        return ss.str();
    }



    static std::string exportVariable(const bool& b)
    {
        return (b == true ? "true" : "false");
    }

    static std::string exportVariable(const std::string& str)
    {
        return '"' + str + '"';
    }


    static std::string exportVariable(const char* str)
    {
        return std::string("\"") + str + '"';
    }


    static std::string exportVariable(const char str)
    {
        return std::string("\"") + str + '"';
    }


    template<typename T1, typename T2>
    static auto exportVariable(const std::pair<T1,T2>& b) -> void
    {
        std::string res;
        res += "{ ";
        res += exportVariable(b.first);
        res += ": ";
        res += exportVariable(b.second);
        res += "}";
    }

    template<typename T>
    static auto exportVariable(const T &list)
      -> typename std::enable_if<
              IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultTraits<T>::value,
             std::string
         >::type
    {
        std::string res;
        auto it = list.begin();
        res += "[ ";
        while (it != list.end()){
            exportVariable(*it);
            if (++it != list.end()){
                res += ", ";
            }
        }
        res += " ]";
    }



    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultTraits<T>::value,
             std::string
         >::type
    {
        std::string res;
        res += "[ ";
        for (decltype (list.size())i = 0; i < list.size(); i++){
            exportVariable(ost, list[i]);
            if (i < list.size() - 1){
                res += ", ";
            }
        }
        res += " ]";
        return res;
    }


    template<typename T>
    static auto exportVariable(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              IsTNLIndexable<T>::value &&
             !IsIndexable<T>::value &&
             !IsIterable<T>::value &&
             !IsExportable<T>::value &&
             !HasDefaultTraits<T>::value,
             std::string
         >::type
    {
        std::string res;
        res += "[ ";
        for (decltype (list.size())i = 0; i < list.size(); i++){
            exportVariable(ost, list[i]);
            if (i <  list.getSize() - 1){
                res += ", ";
            }
        }
        res += " ]";
        return res;
    }



    template<typename T>
    static std::string exportVariable(std::ostream& ost, const std::initializer_list<T> &list)
    {
        std::string res;
        static auto it = list.begin();
        res += "[ ";
        while (it != list.end()){
            exportVariable(ost, *it);
            if (++it != list.end()){
                res += ", ";
            }
        }
        res += " ]";
        return res;
    }


    template<typename T,unsigned int Index = 0, typename Void = void>
    struct PrintClass{
        static std::string print(const T &traitedClass){
            std::string res;
            res += '"' + Traits<T>::traitsType::template getName<Index>() + "\" : ";
            VariableExport::exportVariable(Traits<T>::traitsType::template getReference<Index>()->getValue(traitedClass));
            res += ", ";
            res += PrintClass<T, Index + 1>::print( traitedClass);
            return res;
        }
    };

    template<typename T,unsigned int Index, typename... Types>
    struct PrintClass <Traits<T, Types...>, Index, typename std::enable_if<Index < Traits<T, Types...>::size() - 1>::type>{
        static std::string print(const T &traitedClass){
            std::string res;
            res += '"' + Traits<T, Types...>::template getName<Index>() + "\" : ";
            VariableExport::exportVariable(Traits<T, Types...>::template getReference<Index>()->getValue(traitedClass));
            res += ", ";
            res += PrintClass<Traits<T, Types...>, Index + 1>::print(traitedClass);
            return res;
        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct PrintClass <Traits<T, Types...>, Index, typename std::enable_if<Index == Traits<T, Types...>::size() - 1>::type>{
        static std::string print(const T &traitedClass){
            std::string res;
            res += '"' + Traits<T, Types...>::template getName<Traits<T, Types...>::size() - 1>() + "\" : ";
            res += VariableExport::exportVariable(Traits<T, Types...>::template getReference<Traits<T, Types...>::size() - 1>()->getValue(traitedClass));
            return res;
        }
    };

    template<typename T, unsigned int Index>
    struct PrintClass<T, Index, typename std::enable_if<Index == Traits<T>::traitsType::size() - 1>::type>{
        static std::string print(const T &traitedClass){
            std::string res;
            res += '"' + Traits<T>::traitsType::template getName<Traits<T>::traitsType::size() - 1>() + "\" : ";
            VariableExport::exportVariable(Traits<T>::traitsType::template getReference<Traits<T>::traitsType::size() - 1>()->getValue(traitedClass));
            return res;
        }
    };


    template<typename T>
    static auto exportVariable(const T &traitedClass)
      -> typename std::enable_if<
             HasDefaultIOTraits<T>::value,
             std::string
         >::type
    {
        std::string res;
        res += "{ ";
        res += PrintClass<typename DefaultIOTraits<T>::traitsType>::print(traitedClass);
        res += " }";
        return res;
    }


};

#endif // VARIABLEEXPORT_H
