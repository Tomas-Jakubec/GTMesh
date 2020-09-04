#ifndef PRINTCUSTOM_H
#define PRINTCUSTOM_H
#include "../VariableExport.h"
#include <GTMesh/Traits/CustomTypeTraits.h>

/**
 * @brief Prints an object with defined PrintTo
 * function.
 */
struct PrintCustom {
    static int print(...) {return 0;}

    template <typename T>
    using IsPrintableToStream = Impl::void_t<decltype(PrintTo(std::declval<const T&>(), std::declval<std::ostream&>()))>;


    template <typename T>
    using IsPrintableToStdio = Impl::void_t<decltype(PrintTo(std::declval<const T&>()))>;


    template <typename T>
    static
    IsPrintableToStream<T>
    print(std::ostream& ost, const T& var, ...) {
        PrintTo(var, ost);
    }

    template <typename T>
    static
    IsPrintableToStdio<T>
    print(const T& var, ...) {
        PrintTo(var);
    }
};

#endif // PRINTCUSTOM_H
