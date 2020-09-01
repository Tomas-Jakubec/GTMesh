#ifndef PRINTERS_H
#define PRINTERS_H
#include "../VariableExport.h"


class TraitsPrint {
    int print(...) {}
    template <typename T, std::enable_if_t<HasDefaultIOTraits<T>::value, bool> = true>
    void print(const T& var, std::ostream& ost) {
        VariableExport<VARIABLE_EXPORT_METHOD_OSTREAM>::PrintClass<>
    }
};


#endif // PRINTERS_H
