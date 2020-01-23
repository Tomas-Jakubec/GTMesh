#ifndef TRAITSALGORITHM_H
#define TRAITSALGORITHM_H

#include "../Traits.h"
#include "../CustomTypeTraits.h"
#include <limits>

/*
    operators for arithmetic traits
*/


struct AddProcesor {
    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index < DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    exec(TraitT& op1, const TraitT& op2){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(op1) +=
                DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op2);

        AddProcesor::exec<TraitT, Index + 1>(op1, op2);

    }

    template<typename TraitT, unsigned int Index = 0>
    inline static
    typename std::enable_if<(Index == DefaultArithmeticTraits<TraitT>::size() - 1)>::type
    exec(TraitT& op1,  const TraitT& op2){

        DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getAttr(op1) +=
                DefaultArithmeticTraits<TraitT>::getTraits().template getReference<Index>().getValue(op2);


    }

};

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator+=(TraitT& op1, const TraitT& op2) noexcept {

    AddProcesor::exec(op1, op2);

    return  op1;
}


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(const TraitT& op1, const TraitT& op2) noexcept {
    TraitT res(op1);
    return operator+=(res, op2);
}




template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(const TraitT& op1, TraitT&& op2) noexcept {
    return op2 += op1;
}

/*

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator+(TraitT&& op1, const TraitT& op2) noexcept {
    return op1 += op2;
}
*/


template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator-=(TraitT& op1, const TraitT& op2) noexcept {

    DefaultArithmeticTraits<TraitT>::getTraits().apply(
                [&op1, &op2](unsigned int, const auto& ref, const std::string&) noexcept {
        ref.setValue(op1, ref.getValue(op1) - ref.getValue(op2));
    }
    );
    return  op1;
}

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(const TraitT& op1, const TraitT& op2) noexcept {
    TraitT res(op1);
    return operator-=(res, op2);
}

template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(const TraitT& op1, TraitT&& op2) noexcept {

    return op2 -= op1;
}

/*
template <typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator-(TraitT&& op1, const TraitT& op2) noexcept {

    return op1 -= op2;
}
*/

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const Real& op1, const TraitT& op2) noexcept {
    TraitT res(op2);
    return operator*=(res, op1);
}

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator*(const TraitT& op1, const Real& op2) noexcept {
    TraitT res(op1);
    return operator*=(res, op2);
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

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type&
operator*=(TraitT& op1, const Real& op2) noexcept {

    DefaultArithmeticTraits<TraitT>::getTraits().apply(
                [&op1, &op2](unsigned int, const auto& ref, const std::string&) noexcept {
        ref.setValue(op1, ref.getValue(op1) * op2);
    }
    );
    return  op1;
}


template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator/(const Real& op1, const TraitT& op2) noexcept {

    return (1.0/op1) * op2;
}

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator/(const TraitT& op1, const Real& op2) noexcept {

    return (1.0/op2) * op1;
}


template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator/(const Real& op1, TraitT&& op2) noexcept {

    return (1.0/op1) * op2;
}

template <typename Real, typename TraitT>
typename std::enable_if<HasDefaultArithmeticTraits<TraitT>::value, TraitT>::type
operator/(TraitT&& op1, const Real& op2) noexcept {

    return (1.0/op2) * op1;
}







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

template <typename TraitT>
typename std::enable_if<
    HasDefaultArithmeticTraits<TraitT>::value,
    double //typename TraitCommonType<TraitT>::type
>::type
max(const TraitT& op1) noexcept {
    double res = std::numeric_limits<double>::min();
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
