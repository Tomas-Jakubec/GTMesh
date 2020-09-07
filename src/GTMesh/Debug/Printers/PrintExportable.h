#ifndef PRINTEXPORTABLE_H
#define PRINTEXPORTABLE_H
#include <type_traits>
#include "../VariableExport.h"

struct PrintExportable {
    static int print(...) {return 0;}

    /**
     * Prints the variable with. Valid if the
     * type T has overaloaded the operator<<
     * for ostream.
     */
    template <typename T, std::enable_if_t<IsExportable<T>::value, bool> = true>
    static void print(std::ostream& ost, const T& var, ...) {
        ost << var;
    }

    static void print(std::ostream& ost, const bool& var, ...) {
        ost << (var == true ? "true" : "false");
    }

    static void print(const double& d, ...)
    {
        printf("%g", d);
    }

    static void print(const long double& d, ...)
    {
        printf("%Lg", d);
    }

    static void print(const int& d, ...)
    {
        printf("%d", d);
    }


    static void print(const short& d, ...)
    {
        printf("%hd", d);
    }

    static void print(const unsigned int& d, ...)
    {
        printf("%ud", d);
    }


    static void print(const long long& d, ...)
    {
        printf("%lld", d);
    }

    static void print(const unsigned long long& d, ...)
    {
        printf("%llu", d);
    }

    static void print(const bool& b, ...)
    {
        printf("%s", (b == true ? "true" : "false"));
    }
};

#endif // PRINTEXPORTABLE_H
