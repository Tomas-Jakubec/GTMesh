#ifndef PRINTINDEXABLE_H
#define PRINTINDEXABLE_H
#include "../VariableExport.h"
#include <GTMesh/Traits/CustomTypeTraits.h>
#include <GTMesh/Traits/Interfaces/InterfaceIndexable.h>

struct PrintIndexable {

    template<typename Indexable, std::enable_if_t<Interface::Indexable<Indexable>::value, bool> = true>
    static void print(std::ostream &ost, const Indexable &vec)
    {
        ost << "[ ";
        const auto size = Interface::Indexable<Indexable>::size(vec);
        for (typename std::decay<decltype(size)>::type i = 0; i < size; i++) {
            VariableExport::exportVariable(ost, Interface::Indexable<Indexable>::getElement(vec, i));
            if (i < size - 1) {
                ost << ", ";
            }
        }
        ost << " ]";
    }

    template<typename Indexable,
             typename... TraitsTypes,
             std::enable_if_t<Interface::Indexable<Indexable>::value, bool> = true>
    static void print(std::ostream &ost,
                      const Indexable &vec,
                      const std::tuple<TraitsTypes...> &traitsTuple)
    {
        ost << "[ ";
        const auto size = Interface::Indexable<Indexable>::size(vec);
        for (typename std::decay<decltype(size)>::type i = 0; i < size; i++) {
            VariableExport::exportVariable(ost,
                                           Interface::Indexable<Indexable>::getElement(vec, i),
                                           traitsTuple);
            if (i < size - 1) {
                ost << ", ";
            }
        }
        ost << " ]";
    }

    template<typename Indexable, std::enable_if_t<Interface::Indexable<Indexable>::value, bool> = true>
    static void print(const Indexable &vec)
    {
        printf("[ ");
        const auto size = Interface::Indexable<Indexable>::size(vec);
        for (typename std::decay<decltype(size)>::type i = 0; i < size; i++) {
            VariableExport::exportVariable(Interface::Indexable<Indexable>::getElement(vec, i));
            if (i < size - 1) {
                printf(", ");
            }
        }
        printf(" ]");
    }

    template<typename Indexable,
             typename... TraitsTypes,
             std::enable_if_t<Interface::Indexable<Indexable>::value, bool> = true>
    static void print(const Indexable &vec, const std::tuple<TraitsTypes...> &traitsTuple)
    {
        printf("[ ");
        const auto size = Interface::Indexable<Indexable>::size(vec);
        for (typename std::decay<decltype(size)>::type i = 0; i < size; i++) {
            VariableExport::exportVariable(Interface::Indexable<Indexable>::getElement(vec, i),
                                           traitsTuple);
            if (i < size - 1) {
                printf(", ");
            }
        }
        printf(" ]");
    }
};

#endif // PRINTINDEXABLE_H
