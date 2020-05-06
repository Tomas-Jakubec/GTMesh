#ifndef CONSOLELOGGER_H
#define CONSOLELOGGER_H
#include "VariableExport.h"

#ifdef _WIN32
#ifdef CONSOLE_COLORED_OUTPUT
#include <windows.h>
#endif
#endif
/**
 * @brief The ConsoleLogger class
 */
template <VARIABLE_EXPORT_METHOD method = VARIABLE_EXPORT_METHOD::ostream>
class ConsoleLogger {

public:

    template<typename MSGTYPE, typename... MSGTYPES>
    static void writeMessage(const char* prefix, int line, const char* sourceFile, const MSGTYPE& message, const MSGTYPES&... rest) {
        writeMessage(prefix, line, sourceFile, message);
        writeMessage(prefix, line, sourceFile, rest...);
    }


    template<typename MSGTYPE>
    static void writeMessage(const char* prefix, int line, const char* sourceFile, const MSGTYPE& message) {

#ifdef CONSOLE_COLORED_OUTPUT
#ifdef _WIN32
        HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        std::cerr << prefix << " " << sourceFile << " << " << line << " >> ==> ";
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_INTENSITY);
        VariableExport<>::exportVariable(std::cerr, message);
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
        std::cerr << " <==\n";
#else
        std::cerr << prefix << " " << sourceFile << " << " << line << " >> ==> \033[0;32m";
        VariableExport<>::exportVariable(std::cerr, message);
        std::cerr << "\033[0m <==\n";

#endif
#else
        std::cerr << prefix << " " << sourceFile << " << " << line << " >> ==> ";
        VariableExport<>::exportVariable(std::cerr, message);
        std::cerr << " <==\n";
#endif
    }

    template<typename VAR, typename ... REST>
    static void writeVar(const char* name, VAR value, REST ... rest){

        writeVar(name, value);
        writeVar(rest...);
    }


    template<typename VAR>
    static void writeVar(const char* name, VAR value){

        std::cerr << "variable " << name << " has value: " << value;
    }



    template<typename VAR, typename ... REST>
    static void writeVar(int line, const char* cppFile, const char* name,const VAR& value,const REST& ... rest){

        writeVar(line, cppFile, name, value);
        writeVar(line, cppFile,  rest...);
    }

    template<typename VAR>
    static void writeVar(int line, const char* cppFile, const char* name,const VAR& value){

#ifdef CONSOLE_COLORED_OUTPUT
#ifdef _WIN32
        HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        std::cerr << "== " << cppFile << " << " << line << " >> [[ ";
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
        std::cerr << name;
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
        std::cerr << " ]] ==> ";
        SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
        VariableExport<>::exportVariable(std::cerr, value);
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
        std::cerr << "\n";
#else
        std::cerr << "== " << cppFile << " << " << line << " >> [[ \033[0;33m" << name << "\033[0m ]] ==> \033[0;32m";
        VariableExport<>::exportVariable(std::cerr, value);
        std::cerr << "\033[0m\n";

#endif
#else
        std::cerr << "== " << cppFile << " << " << line << " >> [[ " << name << " ]] ==> ";
        VariableExport<>::exportVariable(std::cerr, value);
        std::cerr << "\n";
#endif
    }

    template<typename VAR>
    static void writeVar(int line, const char* cppFile, const char* name,const std::initializer_list<VAR>& value){

#ifdef CONSOLE_COLORED_OUTPUT
#ifdef _WIN32
        HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        std::cerr << "== " << cppFile << " << " << line << " >> [[ ";
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
        std::cerr << name;
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
        std::cerr << " ]] ==> ";
        SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
        VariableExport<>::exportVariable(std::cerr, value);
        SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
        std::cerr << "\n";
#else
        std::cerr << "== " << cppFile << " << " << line << " >> [[ \033[0;33m" << name << "\033[0m ]] ==> \033[0;32m";
        VariableExport<>::exportVariable(std::cerr, value);
        std::cerr << "\033[0m\n";

#endif
#else
        std::cerr << "== " << cppFile << " << " << line << " >> [[ " << name << " ]] ==> ";
        VariableExport<>::exportVariable(std::cerr, value);
        std::cerr << "\n";
#endif
    }

};




template <>
class ConsoleLogger<VARIABLE_EXPORT_METHOD::stdio> {

public:

    template<typename MSGTYPE, typename... MSGTYPES>
    static void writeMessage(const char* prefix, int line, const char* sourceFile, const MSGTYPE& message, const MSGTYPES&... rest) {
        writeMessage(prefix, line, sourceFile, message);
        writeMessage(prefix, line, sourceFile, rest...);
    }


    template<typename MSGTYPE>
    static void writeMessage(const char* prefix, int line, const char* sourceFile, const MSGTYPE& message) {

        printf("%s %s << %i >> ==> ",
               prefix,
               sourceFile,
               line);

        VariableExport<VARIABLE_EXPORT_METHOD::stdio>::exportVariable(message);

        printf(" <==\n");
    }

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
    static void writeVar(int line, const char* cppFile, const VAR_NAME& name,const VAR& value,const REST& ... rest){

        writeVar(line, cppFile, name, value);
        writeVar(line, cppFile,  rest...);
    }

    template<typename VAR>
    static void writeVar(int line, const char* cppFile, const char* name,const VAR& value){

        printf("== %s << %i >> [[ %s ]] ==> ",
               cppFile,
               line,
               name);

        VariableExport<VARIABLE_EXPORT_METHOD::stdio>::exportVariable(value);
        printf("\n");
    }

    template<typename VAR>
    static void writeVar(int line, const char* cppFile, const char* name,const std::initializer_list<VAR>& value){

        printf("== %s << %i >> [[ %s ]] ==> ",
               cppFile,
               line,
               name);

        VariableExport<VARIABLE_EXPORT_METHOD::stdio>::exportVariable(value);
        printf("\n");
    }

};
#endif // CONSOLELOGGER_H
