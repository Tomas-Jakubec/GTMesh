#ifndef VARIABLEEXPORT_H
#define VARIABLEEXPORT_H
#include <iostream>
#include <fstream>
#include <string>


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
             !Detail::is_exportable<T>::value
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
             !Detail::is_exportable<T>::value
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

};


#endif // VARIABLEEXPORT_H
