/**
 * This file contains the implementation of constexpr for structure.
 * This structure is able to unwind a loop at the compile time.
 */
#ifndef CONSTEXPRFOR_H
#define CONSTEXPRFOR_H
#include <utility>

template <size_t start, size_t stop, size_t step = 1>
class ConstexprFor {
public:
    template<typename Functor,
             size_t Index = start,
             typename ... ArgTypes,
             std::enable_if_t<(start < stop ? Index >= stop : Index <= stop), bool> = true>
    static void exec(ArgTypes&& ...) {}

    template<typename Functor,
             size_t Index = start,
             typename ... ArgTypes,
             std::enable_if_t<(start < stop ? Index < stop : Index > stop), bool> = true>
    static void exec(ArgTypes&& ... args) {
        Functor::template exec<Index>(std::forward<ArgTypes>(args)...);
        exec<Functor, Index + step>(std::forward<ArgTypes>(args)...);
    }

};

template<typename Functor, size_t start, size_t stop, size_t step = 1, typename ... ArgTypes>
void constexprFor(ArgTypes&& ... args) {
    ConstexprFor<start, stop, step>::template exec<Functor>(args...);
}


template<typename Functor, size_t stop, typename ... ArgTypes>
void constexprFor(ArgTypes&& ... args) {
    ConstexprFor<0, stop, 1>::template exec<Functor>(args...);
}

#endif // CONSTEXPRFOR_H
