#ifndef PRINTTUPLE_H
#define PRINTTUPLE_H
#include "../VariableExport.h"
#include <GTMesh/Utils/ConstexprFor.h>
/**
 * @brief Provides functions that prints std::pair and std::tuple similarly to the python.
 */
struct PrintTuple {
    /**
     * @brief Python-like print of dictionary pair of key and value.
     */
    template<typename T1, typename T2>
    static void print(std::ostream& ost, const std::pair<T1,T2>& b)
    {
        ost << "{ ";
        VariableExport::exportVariable(ost, b.first);
        ost << ": ";
        VariableExport::exportVariable(ost, b.second);
        ost << "}";
    }

    template<typename T1, typename T2, typename ... TraitsTypes>
    static void print(std::ostream& ost, const std::pair<T1,T2>& b, const std::tuple<TraitsTypes...>& traitsTuple)
    {
        ost << "{ ";
        VariableExport::exportVariable(ost, b.first, traitsTuple);
        ost << ": ";
        VariableExport::exportVariable(ost, b.second, traitsTuple);
        ost << "}";
    }

    template<typename T1, typename T2>
    static void print(const std::pair<T1,T2>& b)
    {
        printf("{ ");
        VariableExport::exportVariable(b.first);
        printf(": ");
        VariableExport::exportVariable(b.second);
        printf("}");
    }

    template<typename T1, typename T2, typename ... TraitsTypes>
    static void print(const std::pair<T1,T2>& b, const std::tuple<TraitsTypes...>& traitsTuple)
    {
        printf("{ ");
        VariableExport::exportVariable(b.first, traitsTuple);
        printf(": ");
        VariableExport::exportVariable(b.second, traitsTuple);
        printf("}");
    }

    struct PrintTupleExec {
        template < unsigned int Index = 0, typename ... Types>
        static void exec(std::ostream& ost, const std::tuple<Types...>& varTuple) {
            VariableExport::exportVariable(ost, std::get<Index>(varTuple));
            if (Index < sizeof... (Types) - 1) {
                ost << ", ";
            }
        }

        template < unsigned int Index = 0, typename ... Types, typename ... TraitsTypes>
        static void exec(std::ostream& ost, const std::tuple<Types...>& varTuple, const std::tuple<TraitsTypes...>& traitsTuple) {
            VariableExport::exportVariable(ost, std::get<Index>(varTuple), traitsTuple);
            if (Index < sizeof... (Types) - 1) {
                ost << ", ";
            }
        }

        // printf version
        template < unsigned int Index = 0, typename ... Types>
        static void exec(const std::tuple<Types...>& varTuple) {
            VariableExport::exportVariable(std::get<Index>(varTuple));
            if (Index < sizeof... (Types) - 1) {
                printf(", ");
            }
        }

        template < unsigned int Index = 0, typename ... Types, typename ... TraitsTypes>
        static void exec(const std::tuple<Types...>& varTuple, const std::tuple<TraitsTypes...>& traitsTuple) {
            VariableExport::exportVariable(std::get<Index>(varTuple), traitsTuple);
            if (Index < sizeof... (Types) - 1) {
                printf(", ");
            }
        }
    };


    template<typename ... Types>
    static void print(std::ostream& ost, const std::tuple<Types...>& varTuple)
    {
        ost << "[ ";
        constexprFor<PrintTupleExec, sizeof... (Types)>(ost, varTuple);
        ost << "]";
    }


    template<typename ... Types, typename ... TraitsTypes>
    static void print(std::ostream& ost, const std::tuple<Types...>& varTuple, const std::tuple<TraitsTypes...>& traitsTuple)
    {
        ost << "[ ";
        constexprFor<PrintTupleExec, sizeof... (Types)>(ost, varTuple, traitsTuple);
        ost << "]";
    }



    template<typename ... Types>
    static void print(const std::tuple<Types...>& varTuple)
    {
        printf("[ ");
        constexprFor<PrintTupleExec, sizeof... (Types)>(varTuple);
        printf("]");
    }

    template<typename ... Types, typename ... TraitsTypes>
    static void print(const std::tuple<Types...>& varTuple, const std::tuple<TraitsTypes...>& traitsTuple)
    {
        printf("[ ");
        constexprFor<PrintTupleExec, sizeof... (Types)>(varTuple, traitsTuple);
        printf("]");
    }
};

#endif // PRINTTUPLE_H
