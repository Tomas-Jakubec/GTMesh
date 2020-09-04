#ifndef PRINTTUPLE_H
#define PRINTTUPLE_H
#include "../VariableExport.h"

struct PrintTuple {
    static int print(...){return 0;}
    /**
     * @brief Python-like print of dictionary pair of key and value.
     */
    template<typename T1, typename T2>
    static void print(std::ostream& ost, const std::pair<T1,T2>& b, ...)
    {
        ost << "{ ";
        VariableExport<VARIABLE_EXPORT_METHOD_OSTREAM>::exportVariable(ost, b.first);
        ost << ": ";
        VariableExport<VARIABLE_EXPORT_METHOD_OSTREAM>::exportVariable(ost, b.second);
        ost << "}";
    }

    template<typename T1, typename T2>
    static void exportVariable(const std::pair<T1,T2>& b, ...)
    {
        printf("{ ");
        VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(b.first);
        printf(": ");
        VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(b.second);
        printf("}");
    }


    template < unsigned int Index = 0, bool print = false, typename ... Types,  std::enable_if_t<(Index < sizeof... (Types) - 1) && !print, bool > = true >
    static void printTuple(std::ostream& ost, const std::tuple<Types...>& varTuple) {
        printTuple< Index, true >(ost, varTuple);
        printTuple< Index + 1 >(ost, varTuple);
    }


    template < unsigned int Index = 0, bool print = false, typename ... Types,  std::enable_if_t<(Index == sizeof... (Types) - 1) || print, bool > = true >
    static void printTuple(std::ostream& ost, const std::tuple<Types...>& varTuple) {
        VariableExport<VARIABLE_EXPORT_METHOD_OSTREAM>::exportVariable(ost, std::get<Index>(varTuple));
        if (Index < sizeof... (Types) - 1) {
            ost << ", ";
        }
    }

    template<typename ... Types>
    static void print(std::ostream& ost, const std::tuple<Types...>& varTuple, ...)
    {
        ost << "{ ";
        printTuple(ost, varTuple);
        ost << "}";
    }

    // Version for printf
    // These functions are to be performed in CUDA kernels
    template < unsigned int Index = 0, bool print = false, typename ... Types,  std::enable_if_t<(Index < sizeof... (Types) - 1) && !print, bool > = true >
    static void printTuple(const std::tuple<Types...>& varTuple) {
        printTuple< Index, true >(varTuple);
        printTuple< Index + 1 >(varTuple);
    }


    template < unsigned int Index = 0, bool print = false, typename ... Types,  std::enable_if_t<(Index == sizeof... (Types) - 1) || print, bool > = true >
    static void printTuple(const std::tuple<Types...>& varTuple) {
        VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(std::get<Index>(varTuple));
        if (Index < sizeof... (Types) - 1) {
            printf(", ");
        }
    }

    template<typename ... Types>
    static void print(const std::tuple<Types...>& varTuple, ...)
    {
        printf("{ ");
        printTuple(varTuple);
        printf("}");
    }
};

#endif // PRINTTUPLE_H
