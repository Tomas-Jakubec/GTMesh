#ifndef VARIABLEEXPORT_H
#define VARIABLEEXPORT_H
#include <iostream>
#include <fstream>
#include <string>
#include "../Unstructured_mesh/UnstructuredMesh/MeshDataContainer/Traits.h"

namespace Detail {



template <typename T1, typename T2 = void>
struct __is_exportable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_exportable<T1, typename std::enable_if<std::is_class<
        typename std::remove_reference<decltype(std::cerr << std::declval<const T1&>())>::type
        >::value>::type> : public std::integral_constant<bool, true> {

};


template <typename T1>
struct is_exportable : public __is_exportable<T1>
{};


template <typename T1, typename T2 = void>
struct __is_iterable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_iterable<T1, typename std::enable_if<!std::is_same<
        decltype(std::declval<const T1&>().begin()),
        void
     >::value>::type> : public std::integral_constant<bool, true> {

};
template <typename T1>
struct is_iterable : public __is_iterable<T1>
{};


template <typename T1, typename T2 = void>
struct __is_indexable : public std::integral_constant<bool, false> {

};

template <typename T1>
struct __is_indexable<T1, typename std::enable_if<
            !std::is_same<
                decltype (std::declval<const T1&>()[0]), void
            >::value &&
            !std::is_same<
                decltype (std::declval<const T1&>().size()), void
            >::value
    >::type> : public std::integral_constant<bool, true> {

};
template <typename T1>
struct is_indexable : public __is_indexable<T1>
{};



template <typename T1, typename VOID = void>
struct __has_default_traits : public std::integral_constant<bool, false> {
};


template <typename T1>
struct __has_default_traits<
        T1,
        typename std::enable_if<
            Traits<T1>::is_specialized
        >::type
        > : public std::integral_constant<bool, true> {
};

template<typename T>
struct has_default_traits : __has_default_traits<T> {

};
}


struct VariableExport {


    static void _writeWar(std::ostream& ost, ...)
    {
        ost << "\"variable is not exportable\"" << std::endl;
    }


    template<typename T>
    static auto _writeWar(std::ostream& ost, const T& b)
      -> typename std::enable_if<std::is_class<
            typename std::remove_reference<decltype(ost << b)>::type>::value &&
            !std::is_same<T, bool>::value &&
            !std::is_same<T, std::string>::value &&
            !std::is_same<T, const char*>::value &&
            !std::is_same<T, const char>::value
         >::type
    {
        ost << b;
    }



    static void _writeWar(std::ostream& ost, const bool& b)
    {
        ost << (b == true ? "true" : "false");
    }

    static void _writeWar(std::ostream& ost, const std::string& str)
    {
        ost << '"' << str << '"';
    }


    static void _writeWar(std::ostream& ost, const char* str)
    {
        ost << '"' << str << '"';
    }


    static void _writeWar(std::ostream& ost, const char str)
    {
        ost << '"' << str << '"';
    }


    template<typename T1, typename T2>
    static auto _writeWar(std::ostream& ost, const std::pair<T1,T2>& b) -> void
    {
        ost << "{ ";
        _writeWar(ost, b.first);
        ost << ": ";
        _writeWar(ost, b.second);
        ost << "}";
    }

    template<typename T>
    static auto _writeWar(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              Detail::is_iterable<T>::value &&
             !Detail::is_exportable<T>::value &&
             !Detail::has_default_traits<T>::value
         >::type
    {
        auto it = list.begin();
        ost << "[ ";
        while (it != list.end()){
            _writeWar(ost, *it);
            if (++it == list.end()){
                ost << " ]";
            } else {
                ost << ", ";
            }
        }
    }



    template<typename T>
    static auto _writeWar(std::ostream& ost, const T &list)
      -> typename std::enable_if<
              Detail::is_indexable<T>::value &&
             !Detail::is_iterable<T>::value &&
             !Detail::is_exportable<T>::value &&
             !Detail::has_default_traits<T>::value
         >::type
    {
        ost << "[ ";
        for (decltype (list.size())i = 0; i < list.size(); i++){
            _writeWar(ost, list[i]);
            if (i == list.size() - 1){
                ost << " ]";
            } else {
                ost << ", ";
            }
        }
    }






    template<typename T>
    static void _writeWar(std::ostream& ost, const std::initializer_list<T> &list)
    {
        auto it = list.begin();
        ost << "[ ";
        while (it != list.end()){
            _writeWar(ost, *it);
            if (++it == list.end()){
                ost << " ]";
            } else {
                ost << ", ";
            }
        }
    }


    template<typename T,unsigned int Index = 0, typename VOID = void>
    struct PrintClass{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T>::ttype::template getName<Index>() << "\" : ";
            VariableExport::_writeWar(ost, Traits<T>::ttype::template getReference<Index>()->getValue(traitedClass));
            ost << ", ";
            PrintClass<T, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T,unsigned int Index, typename... Types>
    struct PrintClass <Traits<T, Types...>, Index, std::enable_if_t<Index < Traits<T, Types...>::size() - 1>>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T, Types...>::template getName<Index>() << "\" : ";
            VariableExport::_writeWar(ost, Traits<T, Types...>::template getReference<Index>()->getValue(traitedClass));
            ost << ", ";
            PrintClass<Traits<T, Types...>, Index + 1>::print(ost, traitedClass);

        }
    };

    template<typename T,unsigned int Index, typename ... Types>
    struct PrintClass <Traits<T, Types...>, Index, std::enable_if_t<Index == Traits<T, Types...>::size() - 1>>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T, Types...>::template getName<Traits<T, Types...>::size() - 1>() << "\" : ";
            VariableExport::_writeWar(ost, Traits<T, Types...>::template getReference<Traits<T, Types...>::size() - 1>()->getValue(traitedClass));
        }
    };

    template<typename T, unsigned int Index>
    struct PrintClass<T, Index, std::enable_if_t<Index == Traits<T>::ttype::size() - 1>>{
        static void print(std::ostream& ost, const T &traitedClass){
            ost << '"' << Traits<T>::ttype::template getName<Traits<T>::ttype::size() - 1>() << "\" : ";
            VariableExport::_writeWar(ost, Traits<T>::ttype::template getReference<Traits<T>::ttype::size() - 1>()->getValue(traitedClass));
        }
    };


    template<typename T>
    static auto _writeWar(std::ostream& ost, const T &traitedClass)
      -> typename std::enable_if<
             Detail::has_default_traits<T>::value
         >::type
    {
        ost << "{ ";
        PrintClass<typename Traits<T>::ttype>::print(ost, traitedClass);
        ost << " }";
    }


};

#endif // VARIABLEEXPORT_H
