#ifndef MESHFUNCTIONSDEFINE_H
#define MESHFUNCTIONSDEFINE_H

#include <utility>
/**
 * @brief The Order enum declares the order of indexes.
 * The order can be ascend or original.
 */
enum Order{
    ORDER_ASCEND,
    ORDER_ORIGINAL
};

/**
 * @brief The ComputationMethod enum declares the method
 * to be used in mesh properties calculation. The dafault
 * assumes the mesh to have planar elements and
 * star cells with respect to the center.
 * The choice tessellated makes the algorithms more robust.
 * It does not assume planarity and tries to work with the
 * mesh in more generic way.
 */
enum ComputationMethod {
    METHOD_DEFAULT,
    METHOD_TESSELLATED
};


/**
 * @brief The MakeCustomIntegerSequence generates a sequence
 * [StartIndex, StartIndex + Increment, ..., EndIndex]
 */
template <typename Type, Type StartIndex, Type EndIndex, int Increment = 1, Type... Sequence>
struct MakeCustomIntegerSequence : public MakeCustomIntegerSequence<Type, StartIndex + Increment, EndIndex, Increment, Sequence..., StartIndex> {
};

/**
 * @brief The specialization MakeCustomIntegerSequence<Type, EndIndex, EndIndex, Increment, Sequence>
 * terminates the sequence creation
 */
template <typename Type, Type EndIndex, int Increment, Type... Sequence>
struct MakeCustomIntegerSequence<Type, EndIndex, EndIndex, Increment, Sequence...> {
    using type = std::integer_sequence<Type, Sequence..., EndIndex>;
};

/**
 * Alias for MakeCustomIntegerSequence<Type, StartIndex, EndIndex, Increment>::type.
 * See %MakeCustomIntegerSequence.
 */
template<typename Type, Type StartIndex, Type EndIndex, int Increment = 1>
using make_custom_integer_sequence_t = typename MakeCustomIntegerSequence<Type, StartIndex, EndIndex, Increment>::type;


#endif // MESHFUNCTIONSDEFINE_H
