/*
 * A system for parsing function types
 */
#ifndef FUNCTION_TRAITS_H
#define FUNCTION_TRAITS_H
#include <tuple>
template<class F>
struct function_traits;

// function pointer
template<class R, class... Args>
struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)>
{};

template<class R, class... Args>
struct function_traits<R(Args...)>
{
    using return_type = R;

    static constexpr std::size_t arity = sizeof...(Args);

    template <std::size_t N>
    struct argument
    {
        static_assert(N < arity, "error: invalid parameter index.");
        using type = typename std::tuple_element<N,std::tuple<Args...>>::type;
    };

    template <std::size_t N>
    using argument_t = typename argument<N>::type;
};

// member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)> : public function_traits<R(C&,Args...)>
{};

// const member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const> : public function_traits<R(const C&,Args...)>
{};

// member object pointer
template<class C, class R>
struct function_traits<R(C::*)> : public function_traits<R(C&)>
{};

template<typename ...>
using void_t = void;
// body for functor is empty if the operator() ambiguos
// This happens, for example, when the F is the result of std::bind
template <class F, typename = void>
struct function_traits_body{};

template <class F>
struct function_traits_body<F, void_t<decltype(&F::operator())>>
{
    private:
        using call_type = function_traits<decltype(&F::operator())>;
    public:
        using return_type = typename call_type::return_type;

        static constexpr std::size_t arity = call_type::arity - 1;

        template <std::size_t N>
        struct argument
        {
            static_assert(N < arity, "error: invalid parameter index.");
            using type = typename call_type::template argument<N+1>::type;
        };

        template <std::size_t N>
        using argument_t = typename argument<N>::type;
};


// functor
template<class F>
struct function_traits : public function_traits_body<F>
{};

template<class F>
struct function_traits<F&> : public function_traits<F>
{};

template<class F>
struct function_traits<F&&> : public function_traits<F>
{};


#endif // FUNCTION_TRAITS_H
