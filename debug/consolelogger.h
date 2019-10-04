#ifndef CONSOLELOGGER_H
#define CONSOLELOGGER_H
#include <iostream>
#include <fstream>
#include <string>

/* TODO prefer exportable class to iterable
namespace Detail {
constexpr bool is_exportable(...) {
    return false;
}

template <typename T1>
constexpr auto is_exportable(const T1&)
    -> typename std::enable_if<std::is_class<
            typename std::remove_reference<decltype(std::cerr << std::declval<const T1&>())>::type
            >::value
       ,bool>::type
{
    return true;
}
}*/
/**
 * @brief The ConsoleLogger class
 */
class ConsoleLogger {




    static void _writeWar(...)
    {
        std::cerr << "variable is not exportable" << std::endl;
    }


    template<typename T>
    static auto _writeWar(const T& b)
      -> typename std::enable_if<std::is_class<
            typename std::remove_reference<decltype(std::cerr << b)>::type>::value &&
            !std::is_same<T, bool>::value &&
            !std::is_same<T, std::string>::value &&
            !std::is_same<T, const char*>::value
         >::type
    {
        std::cerr << b;
    }



    static void _writeWar(const bool& b)
    {
        std::cerr << (b == true ? "true" : "false");
    }

    static void _writeWar(const std::string& str)
    {
        std::cerr << '"' << str << '"';
    }


    static void _writeWar(const char* str)
    {
        std::cerr << '"' << str << '"';
    }

    template<typename T1, typename T2>
    static auto _writeWar(const std::pair<T1,T2>& b)
    {
        std::cerr << "{ ";
        _writeWar(b.first);
        std::cerr << ", ";
        _writeWar(b.second);
        std::cerr << "}";
    }

    template<typename T>
    static auto _writeWar(const T &list)
      -> typename std::enable_if<
             !std::is_same<
                decltype(std::declval<const T&>().begin()),
                void
             >::value &&
             !std::is_same<T, std::string>::value
         >::type
    {
        auto it = list.begin();
        std::cerr << "[ ";
        while (it != list.end()){
            _writeWar(*it);
            if (++it == list.end()){
                std::cerr << " ]";
            } else {
                std::cerr << ", ";
            }
        }
    }


    template<typename T>
    static void _writeWar(const std::initializer_list<T> &list)
    {
        auto it = list.begin();
        std::cerr << "[ ";
        while (it != list.end()){
            _writeWar(*it);
            if (++it == list.end()){
                std::cerr << " ]";
            } else {
                std::cerr << ", ";
            }
        }
    }


public:



    template<typename VAR_NAME, typename VAR, typename ... REST>
    static void writeVar(VAR_NAME name, VAR value, REST ... rest){

        writeVar(name, value);
        writeVar(rest...);
    }


    template<typename VAR_NAME, typename VAR>
    static void writeVar(VAR_NAME name, VAR value){

        std::cerr << "variable " << name << " has value: " << value;
    }



    template<typename VAR_NAME, typename VAR, typename ... REST>
    static void writeVar(int line, const char* cppFile, VAR_NAME name,const VAR& value, REST ... rest){

        writeVar(line, cppFile, name, value);
        writeVar(line, cppFile,  rest...);
    }

    template<typename VAR_NAME, typename VAR>
    static void writeVar(int line, const char* cppFile, VAR_NAME name,const VAR& value){

#ifdef __linux__
        std::cerr << "In file " << cppFile << " at line " << line << " variable \033[0;33m" << name << "\033[0m has value of \033[0;31m";
        _writeWar(value);
        std::cerr << "\033[0m\n";
#else
        std::cerr << "In file " << cppFile << " at line " << line << " variable " << name << " has value of ";
        _writeWar(value);
        std::cerr << "\n";
#endif
    }

    template<typename VAR_NAME, typename VAR>
    static void writeVar(int line, const char* cppFile, VAR_NAME name,const std::initializer_list<VAR>& value){

#ifdef __linux__
        std::cerr << "In file " << cppFile << " at line " << line << " variable \033[0;33m" << name << "\033[0m has value of \033[0;31m";
        _writeWar(value);
        std::cerr << "\033[0m\n";
#else
        std::cerr << "In file " << cppFile << " at line " << line << " variable " << name << " has value of ";
        _writeWar(value);
        std::cerr << "\n";
#endif
    }

};

#endif // CONSOLELOGGER_H
