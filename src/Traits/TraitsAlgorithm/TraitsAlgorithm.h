#ifndef TRAITSALGORITHM_H
#define TRAITSALGORITHM_H

#include "../Traits.h"
#include "../CustomTypeTraits.h"
#include <limits>
#include <cmath>

/*
    operators for arithmetic traits
*/


/**
 * @brief The TraitsBinaryExpressionProcesor struct applies
 * a binary operator on each member of a traited struct.
 */
template<template<typename, typename> class Operator>
struct TraitsBinaryExpressionProcesor {

/*
    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index < DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& op1, const TraitT& op2){
        Operator<
                typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>
                >::evaluate(
            DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(op1),
            DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op2)
        );
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Index + 1>(op1, op2);

    }

    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index == DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& op1,  const TraitT& op2){

        Operator<
                typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>
                >::evaluate(
            DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(op1),
            DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op2)
        );

    }
*/

    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index < DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& res, const TraitT& op1, const TraitT& op2){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(res) =
                Operator<
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>
                        >::evaluate(
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op1),
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op2)
                );
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Index + 1>(res, op1, op2);

    }

    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index == DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT&res, const TraitT& op1,  const TraitT& op2){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(res) =
                Operator<
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>
                        >::evaluate(
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op1),
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op2)
                );

    }
/*
    template<typename TraitT, typename Real, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index < DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& op1, const Real& op2){
        Operator<
                typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                Real
                >::evaluate(
            DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(op1),
            op2
        );
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Real, Index + 1>(op1, op2);

    }

    template<typename TraitT, typename Real, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index == DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& op1, const Real& op2){
        Operator<
                typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                Real
                >::evaluate(
            DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(op1),
            op2
        );
    }
*/

    template<typename TraitT, typename Real, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index < DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& res, const TraitT& op1, const Real& op2){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(res) =
                Operator<
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                        Real
                        >::evaluate(
                    DefaultArithmeticTraits<TraitT>::getTraits().template getValue<Index>(op1),
                    op2
                );
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Real, Index + 1>(res, op1, op2);

    }

    template<typename TraitT, typename Real, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index == DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT&res, const TraitT& op1,  const Real& op2){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(res) =
                Operator<
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                        Real
                        >::evaluate(
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op1),
                    op2
                );

    }

};



/**
 * @brief The TraitsUnaryExpressionProcesor struct applies
 * a binary operator on each member of a traited struct.
 */
template<template<typename> class Operator>
struct TraitsUnaryExpressionProcesor {


    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index < DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& res, const TraitT& op1){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(res) =
                Operator<
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>
                        >::evaluate(
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op1)
                );
        TraitsUnaryExpressionProcesor<Operator>::evaluate<TraitT, Index + 1>(res, op1);

    }

    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index == DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    evaluate(TraitT& res, const TraitT& op1){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(res) =
                Operator<
                        typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>
                        >::evaluate(
                    DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op1)
                );
    }
};



////
/// Addition operators
///

template< typename T1, typename T2 >
struct BinaryPlus
{
   static auto evaluate( const T1& a, const T2& b ) -> decltype( a + b )
   {
      return a + b;
   }
};

/* This is not necessary since it is equivalent to a = a + b
template< typename T1, typename T2 >
struct PlusEqual
{
   static auto evaluate( T1& a, const T2& b ) -> decltype( a += b )
   {
      return a += b;
   }
};
*/

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator+=(TraitT& op1, const TraitT& op2) noexcept {

    TraitsBinaryExpressionProcesor<BinaryPlus>::evaluate(op1, op1, op2);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(const TraitT& op1, const TraitT& op2) noexcept {

    TraitT res;
    TraitsBinaryExpressionProcesor<BinaryPlus>::evaluate(res, op1, op2);
    return res;
}




template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(const TraitT& op1, TraitT&& op2) noexcept {

    return op2 += op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(TraitT&& op1, const TraitT& op2) noexcept {
    return op1 += op2;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(TraitT&& op1, TraitT&& op2) noexcept {
    return op2 += op1;
}


////
/// Subtraction operators
///

template< typename T1, typename T2 >
struct BinaryMinus
{
   static auto evaluate( const T1& a, const T2& b ) -> decltype( a - b )
   {
      return a - b;
   }
};


/*This is not necessary since it is equivalent to a = a - b
template< typename T1, typename T2 >
struct MinusEqual
{
   static auto evaluate( T1& a, const T2& b ) -> decltype( a -= b )
   {
      return a -= b;
   }
};
*/

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator-=(TraitT& op1, const TraitT& op2) noexcept {

    TraitsBinaryExpressionProcesor<BinaryMinus>::evaluate(op1, op1, op2);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(const TraitT& op1, const TraitT& op2) noexcept {

    TraitT res;
    TraitsBinaryExpressionProcesor<BinaryMinus>::evaluate(res, op1, op2);
    return res;
}




template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(const TraitT& op1, TraitT&& op2) noexcept {

    return op2 -= op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(TraitT&& op1, const TraitT& op2) noexcept {
    return op1 -= op2;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(TraitT&& op1, TraitT&& op2) noexcept {
    return op2 -= op1;
}





////
/// Multiplication operators
///

template< typename T1, typename T2 >
struct BinaryStar
{
   static auto evaluate( const T1& a, const T2& b ) -> decltype( a * b )
   {
      return a * b;
   }
};

/* This is not necessary since it is equivalent to a = a * b
template< typename T1, typename T2 >
struct StarEqual
{
   static auto evaluate( T1& a, const T2& b ) -> decltype( a *= b )
   {
      return a *= b;
   }
};
*/


///
/// Multiplication of traited class with a real number
///
template <typename TraitT, typename Real>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator*=(TraitT& op1, const Real& op2) noexcept {

    TraitsBinaryExpressionProcesor<BinaryStar>::evaluate(op1, op1, op2);

    return  op1;
}


template <typename TraitT, typename Real>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const TraitT& op1, const Real& op2) noexcept {

    TraitT res;
    TraitsBinaryExpressionProcesor<BinaryStar>::evaluate(res, op1, op2);
    return res;
}


template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const Real& op1, const TraitT& op2) noexcept {
    return operator*(op2, op1);
}



template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const Real& op1, TraitT&& op2) noexcept {
    return operator*=(op2, op1);
}

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(TraitT&& op1, const Real& op2) noexcept {
    return operator*=(op1, op2);
}




///
/// Multiplication of two traited classes (by elements)
///

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator*=(TraitT& op1, const TraitT& op2) noexcept {

    TraitsBinaryExpressionProcesor<BinaryStar>::evaluate(op1, op1, op2);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const TraitT& op1, const TraitT& op2) noexcept {

    TraitT res;
    TraitsBinaryExpressionProcesor<BinaryStar>::evaluate(res, op1, op2);
    return res;
}

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const TraitT& op1, TraitT&& op2) noexcept {
    return operator*=(op2, op1);
}

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(TraitT&& op1, const TraitT& op2) noexcept {
    return operator*=(op1, op2);
}


///
/// Division of traited class with a real number
///

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator/(const TraitT& op1, const Real& op2) noexcept {

    return (1.0/op2) * op1;
}


template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator/(TraitT&& op1, const Real& op2) noexcept {

    return (1.0/op2) * op1;
}



////
/// Power
///

template< typename T1, typename T2 >
struct Pow
{
   static auto evaluate( const T1& a, const T2& b ) -> decltype( pow(a,b) )
   {
      return pow(a,b);
   }
};




///
/// Power of traited class to a real number
///
template <typename TraitT, typename Real>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
pow(TraitT&& op1, const Real& op2) noexcept {

    TraitsBinaryExpressionProcesor<Pow>::evaluate(op1, op1, op2);

    return  op1;
}


template <typename TraitT, typename Real>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
pow(const TraitT& op1, const Real& op2) noexcept {

    TraitT res;
    TraitsBinaryExpressionProcesor<Pow>::evaluate(res, op1, op2);
    return res;
}


template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
pow(const Real& op1, const TraitT& op2) noexcept {
    return pow(op2, op1);
}



template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
pow(const Real& op1, TraitT&& op2) noexcept {
    return pow(op2, op1);
}


/*
template <typename T, typename Real>
typename std::enable_if<
    !std::is_class<T>::value,
    T
>::type
pow(const T& arg, const Real& p) noexcept {

    return  std::pow(arg, p);
}
*/
template <typename T, typename Real>
typename std::enable_if<
    IsIndexable<T>::value,
    T
>::type
pow(const T& array, const Real& p) noexcept {

    T resArray(array);

    for (decltype (resArray.size()) index = 0; index < resArray.size(); index++){
        resArray[index] = pow(resArray[index], p);

    }
    return  resArray;
}

template <typename T,typename Real>
typename std::enable_if<
    IsTNLIndexable<T>::value,
    T
>::type
pow(const T& array, const Real& p) noexcept {

    T resArray(array);

    for (decltype (resArray.getSize()) index = 0; index < resArray.getSize(); index++){
        resArray[index] = pow(resArray[index], p);

    }
    return  resArray;
}


////
/// Exp
///

template< typename T1 >
struct Exp
{
   static auto evaluate( const T1& a ) -> decltype( std::exp(a) )
   {
      return std::exp(a);
   }
};




///
/// E to each traited class member
///
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
exp(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Exp>::evaluate(op1, op1);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
exp(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Exp>::evaluate(res, op1);
    return res;
}





////
/// Log
///

template< typename T1 >
struct Log
{
   static auto evaluate( const T1& a ) -> decltype( std::log(a) )
   {
      return std::log(a);
   }
};




///
/// Logarithm of each traited class member
///
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
log(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Log>::evaluate(op1, op1);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
log(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Log>::evaluate(res, op1);
    return res;
}



////
/// Log2
///

template< typename T1 >
struct Log2
{
   static auto evaluate( const T1& a ) -> decltype( std::log2(a) )
   {
      return std::log2(a);
   }
};




///
/// Logarithm of each traited class member
///
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
log2(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Log2>::evaluate(op1, op1);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
log2(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Log2>::evaluate(res, op1);
    return res;
}




////
/// Log10
///

template< typename T1 >
struct Log10
{
   static auto evaluate( const T1& a ) -> decltype( std::log10(a) )
   {
      return std::log10(a);
   }
};




///
/// Logarithm of each traited class member
///
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
log10(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Log10>::evaluate(op1, op1);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
log10(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Log10>::evaluate(res, op1);
    return res;
}



////
/// Sqrt
///

template< typename T1 >
struct Sqrt
{
   static auto evaluate( const T1& a ) -> decltype( std::sqrt(a) )
   {
      return std::sqrt(a);
   }
};




///
/// Sqrt of each traited class member
///
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
sqrt(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Sqrt>::evaluate(op1, op1);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
sqrt(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Sqrt>::evaluate(res, op1);
    return res;
}


////
/// Cbrt
///

template< typename T1 >
struct Cbrt
{
   static auto evaluate( const T1& a ) -> decltype( std::cbrt(a) )
   {
      return std::cbrt(a);
   }
};




///
/// Cbrt of each traited class member
///
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
cbrt(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Cbrt>::evaluate(op1, op1);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
cbrt(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Cbrt>::evaluate(res, op1);
    return res;
}

///
/// Unary minus
///

template< typename T1 >
struct UnaryMinus
{
   static auto evaluate( const T1& a ) -> decltype( -1 * a )
   {
      return -a;
   }
};

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<UnaryMinus>::evaluate(res, op1);
    return res;
}

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator-(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<UnaryMinus>::evaluate(op1, op1);

    return op1;
}


///
/// Absolute value
///

template< typename T1 >
struct Abs
{

    static auto evaluate( const T1& a ) -> T1 {
        return abs(a);
    }

};




template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
abs(const TraitT& op1) noexcept {

    TraitT res;
    TraitsUnaryExpressionProcesor<Abs>::evaluate(res, op1);
    return res;
}

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
abs(TraitT&& op1) noexcept {

    TraitsUnaryExpressionProcesor<Abs>::evaluate(op1, op1);

    return op1;
}



template <typename T>
typename std::enable_if<
    !std::is_class<T>::value,
    T
>::type
abs(const T& arg) noexcept {

    return  std::abs(arg);
}

template <typename T>
typename std::enable_if<
    IsIndexable<T>::value,
    T
>::type
abs(const T& array) noexcept {

    T resArray(array);

    for (decltype (resArray.size()) index = 0; index < resArray.size(); index++){
        resArray[index] = abs(resArray[index]);

    }
    return  resArray;
}
/*
template <typename T>
typename std::enable_if<
    IsTNLIndexable<T>::value,
    T
>::type
abs(const T& array) noexcept {

    T resArray(array);

    for (decltype (resArray.getSize()) index = 0; index < resArray.getSize(); index++){
        resArray[index] = abs(resArray[index]);

    }
    return  resArray;
}
*/







template<typename Trait, typename Void = void>
class TraitCommonType
{
public:
    using type = void;
};

/*
template<typename Trait>
class TraitCommonType
        <
        Trait,
        typename std::enable_if<
        HasDefaultArithmeticTraits<Trait>::value &&
        DefaultArithmeticTraits<Trait>::traitType::size() == 1
        >::type
        >
{
public:
    using type = typename DefaultArithmeticTraits<Trait>::template type<0>;
};*/

template<typename Trait>
class TraitCommonType<Trait,
        typename std::enable_if<HasDefaultArithmeticTraits<Trait>::value /*&&(DefaultArithmeticTraits<Trait>::traitType::size() > 1)*/>::type
>
{
    using traitType = DefaultArithmeticTraits<Trait>;


    template<unsigned int Index = 0, typename Void = void, typename Dummy = void>
    struct declComonType{
        static_assert (sizeof (Void) != 0, "Tudy vlastne chci projit");
        using type = typename std::common_type<typename declComonType<Index + 1>::type, typename traitType::template type<Index>>::type;
    };
    /*template<unsigned int Index, typename Dummy>
    struct declComonType<Index,
            typename std::enable_if<HasDefaultArithmeticTraits<typename traitType::template type<Index>>::value>::type,
            Dummy>{
        using type = typename std::common_type<typename TraitCommonType<typename traitType::template type<Index>>::type, typename declComonType<Index + 1>::type>::type;
    };*/
    template<typename Dummy>
    struct declComonType<traitType::traitType::size() - 1,
            typename std::enable_if<!HasDefaultArithmeticTraits<typename traitType::template type<traitType::traitType::size() - 1>>::value>::type
            , Dummy>{

        static_assert (sizeof(Dummy) != 0, "Tudy vlastne chci projit");
        using type = typename traitType::template type<traitType::traitType::size() - 1>;
    };
/*
    template<typename Dummy>
    struct declComonType<traitType::traitType::size() - 1,
            typename std::enable_if<HasDefaultArithmeticTraits<typename traitType::template type<traitType::traitType::size() - 1>>::value>::type,
            Dummy>{
        using type = typename TraitCommonType<typename traitType::template type<traitType::traitType::size() - 1>>::type;
    };
*/
public:
    using type = typename declComonType<>::type;
};




template <typename T>
typename std::enable_if<
    !std::is_class<T>::value,
    T
>::type
max(const T& arg) noexcept {

    return  arg;
}

template <typename T>
typename std::enable_if<
    IsIndexable<T>::value,
    double
>::type
max(const T& array) noexcept {

    double res = std::numeric_limits<double>::min();
    for (decltype (array.size()) index = 0; index < array.size(); index++){
        double m = max(array[index]);
        if (res < m){
            res = m;
        }

    }
    return  res;
}

template <typename T>
typename std::enable_if<
    IsTNLIndexable<T>::value,
    double
>::type
max(const T& array) noexcept {

    double res = std::numeric_limits<double>::min();
    for (decltype (array.getSize()) index = 0; index < array.getSize(); index++){
        double m = max(array[index]);
        if (res < m){
            res = m;
        }

    }
    return  res;
}

template <typename TraitT>
typename std::enable_if<
    HasDefaultArithmeticTraits<TraitT>::value,
    double //typename TraitCommonType<TraitT>::type
>::type
max(const TraitT& op1) noexcept {
    double res = -std::numeric_limits<double>::max();
    //typename TraitCommonType<TraitT>::type res = std::numeric_limits<typename TraitCommonType<TraitT>::type>::min();

    DefaultArithmeticTraits<TraitT>::getTraits().apply(
                [&res, &op1](unsigned int, const auto& ref, const std::string&) noexcept {
         auto m = max(ref.getValue(op1));
         if (res < m) {
             res = m;
         }
    }
    );
    return  res;
}



template <typename T>
typename std::enable_if<
    !std::is_class<T>::value,
    T
>::type
min(const T& arg) noexcept {

    return  arg;
}

template <typename T>
typename std::enable_if<
    IsIndexable<T>::value,
    double
>::type
min(const T& array) noexcept {

    double res = std::numeric_limits<double>::max();
    for (decltype (array.size()) index = 0; index < array.size(); index++){
        double m = min(array[index]);
        if (res > m){
            res = m;
        }

    }
    return  res;
}

template <typename TraitT>
typename std::enable_if<
    HasDefaultArithmeticTraits<TraitT>::value,
    double//typename TraitCommonType<TraitT>::type
>::type
min(const TraitT& op1) noexcept {
    double res = std::numeric_limits<double>::max();
    //typename TraitCommonType<TraitT>::type res = std::numeric_limits<typename TraitCommonType<TraitT>::type>::max();

    DefaultArithmeticTraits<TraitT>::getTraits().apply(
                [&res, &op1](unsigned int, const auto& ref, const std::string&) noexcept {
         auto m = min(ref.getValue(op1));
         if (res > m) {
             res = m;
         }
    }
    );
    return  res;
}







#endif // TRAITSALGORITHM_H


