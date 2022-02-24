#ifndef TRAITSALGORITHM_H
#define TRAITSALGORITHM_H

#include "../Traits.h"
#include "../CustomTypeTraits.h"
#include <limits>
#include <cmath>
#include <algorithm>

/**
 * The TraitsAlhorithm is an extension of Traits. It provides automatically generated
 * operators and other mathematic operations,
 * (e.g. addition, subtraction, multiplication, maximum) for the classes with DefaultArithmeticTraits defined.
 * The operations are implemented element wise. (The classes are treated as small static vectors of data)
 *
 * So far, the algorithms supports traits with direct refferences only.
*/

/**
 * @brief The TraitsBinaryExpressionProcesor struct applies
 * a binary operator on each member of a traited struct.
 */
template<template<typename, typename> class Operator>
struct TraitsBinaryExpressionProcesor
{
    template<typename TraitT,
             unsigned int Index = 0,
             bool ApplyOperation = Index == DefaultArithmeticTraits<TraitT>::size() - 1>
    inline static typename std::enable_if<!ApplyOperation>::type evaluate(TraitT &res,
                                                                          const TraitT &op1,
                                                                          const TraitT &op2)
    {
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Index, true>(res, op1, op2);
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Index + 1>(res, op1, op2);
    }

    template<typename TraitT,
             unsigned int Index = 0,
             bool ApplyOperation = Index == DefaultArithmeticTraits<TraitT>::size() - 1>
    inline static typename std::enable_if<ApplyOperation>::type evaluate(TraitT &res,
                                                                         const TraitT &op1,
                                                                         const TraitT &op2)
    {
        DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(res)
            = Operator<typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                       typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>>::
                evaluate(DefaultArithmeticTraits<TraitT>::getTraits().template getValue<Index>(op1),
                         DefaultArithmeticTraits<TraitT>::getTraits().template getValue<Index>(op2));
    }

    template<typename TraitT,
             typename Real,
             unsigned int Index = 0,
             bool ApplyOperation = Index == DefaultArithmeticTraits<TraitT>::size() - 1>
    inline static typename std::enable_if<!ApplyOperation>::type evaluate(TraitT &res,
                                                                          const TraitT &op1,
                                                                          const Real &op2)
    {
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Real, Index, true>(res, op1, op2);
        TraitsBinaryExpressionProcesor<Operator>::evaluate<TraitT, Real, Index + 1>(res, op1, op2);
    }

    template<typename TraitT,
             typename Real,
             unsigned int Index = 0,
             bool ApplyOperation = Index == DefaultArithmeticTraits<TraitT>::size() - 1>
    inline static typename std::enable_if<ApplyOperation>::type evaluate(TraitT &res,
                                                                         const TraitT &op1,
                                                                         const Real &op2)
    {
        DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(res)
            = Operator<typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>,
                       Real>::evaluate(DefaultArithmeticTraits<TraitT>::getTraits()
                                           .template getValue<Index>(op1),
                                       op2);
    }
};

/**
 * @brief The TraitsUnaryExpressionProcesor struct applies
 * a binary operator on each member of a traited struct.
 */
template<template<typename> class Operator>
struct TraitsUnaryExpressionProcesor
{
    template<typename TraitT,
             unsigned int Index = 0,
             bool ApplyOperation = Index == DefaultArithmeticTraits<TraitT>::size() - 1>
    inline static typename std::enable_if<!ApplyOperation>::type evaluate(TraitT &res,
                                                                          const TraitT &op1)
    {
        TraitsUnaryExpressionProcesor<Operator>::evaluate<TraitT, Index, true>(res, op1);
        TraitsUnaryExpressionProcesor<Operator>::evaluate<TraitT, Index + 1>(res, op1);
    }

    template<typename TraitT,
             unsigned int Index = 0,
             bool ApplyOperation = Index == DefaultArithmeticTraits<TraitT>::size() - 1>
    inline static typename std::enable_if<ApplyOperation>::type evaluate(TraitT &res,
                                                                         const TraitT &op1)
    {
        DefaultArithmeticTraits<TraitT>::getTraits().template getAttr<Index>(res)
            = Operator<typename DefaultArithmeticTraits<TraitT>::traitsType::template type<Index>>::
                evaluate(DefaultArithmeticTraits<TraitT>::getTraits().template getValue<Index>(op1));
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

    TraitsBinaryExpressionProcesor<BinaryMinus>::evaluate(op2, op1, op2);

}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(TraitT&& op1, const TraitT& op2) noexcept {
    return op1 -= op2;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(TraitT&& op1, TraitT&& op2) noexcept {
    return op1 -= op2;
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
/// Multiplication of two traited classes (elementwise)
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
   static auto evaluate( const T1& a ) -> decltype( sqrt(a) )
   {
        using std::sqrt;
        return sqrt(a);
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
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
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

    return  abs(arg);
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







template<typename Trait, typename Operation, typename = void>
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

template<typename Trait, typename Operation>
class TraitCommonType<Trait, Operation,
        typename std::enable_if<HasDefaultArithmeticTraits<Trait>::value /*&&(DefaultArithmeticTraits<Trait>::traitType::size() > 1)*/>::type
>
{
    using traitType = typename DefaultArithmeticTraits<Trait>::traitsType;


    template<unsigned int Index = 0, typename Void = void, typename Dummy = void>
    struct declComonType{
        static_assert (sizeof (Void) != 0, "Tudy vlastne chci projit");
        using type = typename std::common_type<typename declComonType<Index + 1>::type, decltype (Operation(std::declval<typename traitType::template type<Index>>()))>::type;
    };

    template<unsigned int Index, typename Dummy>
    struct declComonType<Index,
            typename std::enable_if<Index == traitType::size() - 1>::type
            , Dummy>{

        static_assert (sizeof(Dummy) != 0, "Tudy vlastne chci projit");
        using type = decltype (Operation(std::declval<typename traitType::template type<Index>>()));
    };

public:
    using type = typename declComonType<>::type;
};

template<template<typename, typename> class Operator>
struct TraitsAggregationProcesor {
private:
    template <typename TraitType>
    static constexpr
        std::enable_if_t<HasDefaultArithmeticTraits<TraitType>::value, unsigned int>
        lastAttributeIndex(){
        return DefaultArithmeticTraits<TraitType>::size() - 1;
    }

    template <typename TraitType>
    static constexpr
        std::enable_if_t<!HasDefaultArithmeticTraits<TraitType>::value, unsigned int>
        lastAttributeIndex(...){
        return 0;
    }

public:
    template<typename TraitT,
             unsigned int Index = lastAttributeIndex<TraitT>(),
             typename std::enable_if<Index == 0 && HasDefaultArithmeticTraits<TraitT>::value, bool>::type = true>
    inline static auto evaluate(const TraitT &op1)
    {
        return evaluate(DefaultArithmeticTraits<TraitT>::getTraits().template getValue<Index>(op1));
    }

    template<typename TraitT,
             unsigned int Index = lastAttributeIndex<TraitT>(),
             typename std::enable_if<(Index > 0) && (Index <= lastAttributeIndex<TraitT>())
                                         && HasDefaultArithmeticTraits<TraitT>::value,
                                     bool>::type = true>
    inline static auto evaluate(const TraitT &op1)
    {
        return Operator<
            decltype(evaluate<TraitT, Index - 1>(op1)),
            decltype(evaluate(DefaultArithmeticTraits<TraitT>::getTraits().template getValue<Index>(
                op1)))>::evaluate((evaluate<TraitT, Index - 1>(op1)),
                                  (evaluate(DefaultArithmeticTraits<TraitT>::getTraits()
                                                .template getValue<Index>(op1))));
    }

    template <typename T, typename std::enable_if< !std::is_class<T>::value, bool >::type = true>
    inline static
        auto
        evaluate(const T& arg) noexcept {
        return  arg;
    }

    template <typename T, typename std::enable_if< IsIndexable<T>::value, bool >::type = true>
    inline static
        auto
        evaluate(const T& array) noexcept {
        if (array.size() > 0){
            using resType = decltype (evaluate(array[0]));
            auto res = evaluate(array[0]);
            for (decltype (array.size()) index = 1; index < array.size(); index++){
                res = Operator<resType, resType>::evaluate(res, evaluate(array[index]));
            }

            return  res;
        }
        return 0.0;
    }


};

template <typename T1, typename T2>
struct Max{
    static auto evaluate(const T1& _1, const T2& _2){
        return std::max<std::common_type_t<T1, T2>>(_1, _2);
    }
};

template <typename T>
auto max(const T& val){
    return TraitsAggregationProcesor<Max>::evaluate(val);
}


template <typename T1, typename T2>
struct Min{
    static auto evaluate(const T1& _1, const T2& _2){
        return std::min<std::common_type_t<T1, T2>>(_1, _2);
    }
};

template <typename T>
auto min(const T& val){
    return TraitsAggregationProcesor<Min>::evaluate(val);
}


template <typename T>
auto sum(const T& val){
    return TraitsAggregationProcesor<BinaryPlus>::evaluate(val);
}






#endif // TRAITSALGORITHM_H



