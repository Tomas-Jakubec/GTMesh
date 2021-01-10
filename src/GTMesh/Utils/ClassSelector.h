#ifndef CLASSSELECTOR_H
#define CLASSSELECTOR_H
#include <type_traits>


template<typename VarType, template<typename, typename> class Condition, typename Class, typename... Printers>
struct ClassSelector : public std::conditional_t<Condition<VarType, Class>::value,
                                                 ClassSelector<VarType, Condition, Class>,
                                                 ClassSelector<VarType, Condition, Printers...>>
{};

template<typename VarType,template<typename, typename> class Condition, typename Class>
struct ClassSelector<VarType, Condition, Class>
{
    using SelectedClass = Class;
};

#endif // CLASSSELECTOR_H
