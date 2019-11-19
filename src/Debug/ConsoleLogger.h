#ifndef CONSOLELOGGER_H
#define CONSOLELOGGER_H
#include "VariableExport.h"

// TODO prefer exportable class to iterable
/**
 * @brief The ConsoleLogger class
 */
class ConsoleLogger {

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
        VariableExport::_writeWar(std::cerr, value);
        std::cerr << "\033[0m\n";
#else
        std::cerr << "In file " << cppFile << " at line " << line << " variable " << name << " has value of ";
        VariableExport::exportVariable(std::cerr, value);
        std::cerr << "\n";
#endif
    }

    template<typename VAR_NAME, typename VAR>
    static void writeVar(int line, const char* cppFile, VAR_NAME name,const std::initializer_list<VAR>& value){

#ifdef __linux__
        std::cerr << "In file " << cppFile << " at line " << line << " variable \033[0;33m" << name << "\033[0m has value of \033[0;31m";
        VariableExport::_writeWar(std::cerr, value);
        std::cerr << "\033[0m\n";
#else
        std::cerr << "In file " << cppFile << " at line " << line << " variable " << name << " has value of ";
        VariableExport::exportVariable(std::cerr, value);
        std::cerr << "\n";
#endif
    }

};

#endif // CONSOLELOGGER_H
