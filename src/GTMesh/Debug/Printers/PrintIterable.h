#ifndef PRINTITERABLE_H
#define PRINTITERABLE_H
#include "../VariableExport.h"
#include <GTMesh/Traits/CustomTypeTraits.h>

struct PrintIterable {

    template<typename T>
    using isIterable = std::enable_if_t < IsIterable<T>::value, bool >;

    template< typename Iterable,
              isIterable<Iterable> = true >
    static void print(std::ostream& ost, const Iterable& list) {
        auto it = list.begin();
        ost << "[ ";
        while (it != list.end()){
            VariableExport<>::exportVariable(ost, *it);
            if (++it != list.end()){
                ost << ", ";
            }
        }
        ost << " ]";
    }

    template< typename Iterable,
              typename ... TraitsTypes,
              isIterable<Iterable> = true >
    static void print(std::ostream& ost, const Iterable& list, const std::tuple<TraitsTypes...>& traitsTuple) {
        auto it = list.begin();
        ost << "[ ";
        while (it != list.end()){
            VariableExport<>::exportVariable(ost, *it, traitsTuple);
            if (++it != list.end()){
                ost << ", ";
            }
        }
        ost << " ]";
    }

    template< typename Iterable,
              isIterable<Iterable> = true >
    static void print(const Iterable &list)
    {
        auto it = list.begin();
        printf("[ ");
        while (it != list.end()){
            VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(*it);
            if (++it != list.end()){
                printf(", ");
            }
        }
        printf(" ]");
    }

    template< typename Iterable,
              typename ... TraitsTypes,
              isIterable<Iterable> = true >
    static void print(const Iterable &list, const std::tuple<TraitsTypes...>& traitsTuple)
    {
        auto it = list.begin();
        printf("[ ");
        while (it != list.end()){
            VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(*it, traitsTuple);
            if (++it != list.end()){
                printf(", ");
            }
        }
        printf(" ]");
    }
};

#endif // PRINTITERABLE_H
