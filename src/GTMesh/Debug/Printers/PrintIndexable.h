#ifndef PRINTINDEXABLE_H
#define PRINTINDEXABLE_H
#include "../VariableExport.h"
#include <GTMesh/Traits/CustomTypeTraits.h>

struct PrintIndexable {

    template <typename Indexable, std::enable_if_t<IsIndexable<Indexable>::value, bool> = true>
    static auto size(const Indexable& vec) {
        return vec.size();
    }

    template <typename Indexable, std::enable_if_t<IsTNLIndexable<Indexable>::value, bool> = true>
    static auto size(const Indexable& vec) {
        return vec.getSize();
    }

    template< typename Indexable, std::enable_if_t<IsIndexable<Indexable>::value || IsTNLIndexable<Indexable>::value, bool> = true >
    static void print(std::ostream& ost, const Indexable& vec) {
        ost << "[ ";
        for (decltype (size(vec))i = 0; i < size(vec); i++){
            VariableExport<>::exportVariable(ost, vec[i]);
            if (i <  size(vec) - 1){
                ost << ", ";
            }
        }
        ost << " ]";
    }


    template< typename Indexable, typename ... TraitsTypes, std::enable_if_t<IsIndexable<Indexable>::value || IsTNLIndexable<Indexable>::value, bool> = true >
    static void print(std::ostream& ost, const Indexable& vec, const std::tuple<TraitsTypes...>& traitsTuple) {
        ost << "[ ";
        for (decltype (size(vec))i = 0; i < size(vec); i++){
            VariableExport<>::exportVariable(ost, vec[i], traitsTuple);
            if (i <  size(vec) - 1){
                ost << ", ";
            }
        }
        ost << " ]";
    }


    template< typename Indexable, std::enable_if_t<IsIndexable<Indexable>::value || IsTNLIndexable<Indexable>::value, bool> = true >
    static void print(const Indexable& vec) {
        printf("[ ");
        for (decltype (size(vec))i = 0; i < size(vec); i++){
            VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(vec[i]);
            if (i < size(vec) - 1){
                printf(", ");
            }
        }
        printf(" ]");
    }


    template< typename Indexable, typename ... TraitsTypes, std::enable_if_t<IsIndexable<Indexable>::value || IsTNLIndexable<Indexable>::value, bool> = true >
    static void print(const Indexable& vec, const std::tuple<TraitsTypes...>& traitsTuple) {
        printf("[ ");
        for (decltype (size(vec))i = 0; i < size(vec); i++){
            VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(vec[i], traitsTuple);
            if (i < size(vec) - 1){
                printf(", ");
            }
        }
        printf(" ]");
    }
};

#endif // PRINTINDEXABLE_H
