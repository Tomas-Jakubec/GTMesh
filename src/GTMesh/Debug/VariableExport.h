#ifndef VARIABLEEXPORT_H
#define VARIABLEEXPORT_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../Traits/Traits.h"
#include "../Traits/CustomTypeTraits.h"
#include "../Traits/TraitsBind/TraitsBind.h"
#include "../Traits/CustomTypeTraits.h"




struct VariableExport {


    template <typename T>
    static void exportVariable(std::ostream& ost, const T& var);


    template <typename T, typename ... TraitsTypes>
    static void exportVariable(std::ostream& ost, const T& var, const std::tuple<TraitsTypes...>&);


    template <typename T>
    static void exportVariable(const T& var);


    template <typename T, typename ... TraitsTypes>
    static void exportVariable(const T& var, const std::tuple<TraitsTypes...>&);

};


/*
 * Printers
 */
#include <GTMesh/Utils/ClassSelector.h>
#include "Printers/PrintAny.h"
#include "Printers/PrintCustom.h"
#include "Printers/PrintExportable.h"
#include "Printers/PrintText.h"
#include "Printers/PrintIterable.h"
#include "Printers/PrintIndexable.h"
#include "Printers/PrintTuple.h"
#include "Printers/PrintTraitedClass.h"




namespace Impl {
template<typename VarType, typename Printer, typename = void>
struct IsPrintableBy : public std::false_type
{};

template<typename VarType, typename Printer>
struct IsPrintableBy<VarType,
                     Printer,
                     decltype(Printer::print(std::declval<std::ostream &>(),
                                             std::declval<const VarType &>()))>
    : public std::true_type
{};

template<typename VarType, typename TraitsTuple, typename Printer, typename = void>
struct IsPrintableByWithTuple : public std::false_type
{};

template<typename VarType, typename TraitsTuple, typename Printer>
struct IsPrintableByWithTuple<VarType,
                              TraitsTuple,
                              Printer,
                              decltype(Printer::print(std::declval<std::ostream &>(),
                                                      std::declval<const VarType &>(),
                                                      std::declval<const TraitsTuple &>()))>
    : public std::true_type
{};

} // Impl namespace

/// @brief This class detects whether the given Printer class has a suitable print function
/// to export the variable of VarType.
template<typename VarType, typename Printer>
struct IsPrintableBy: public Impl::IsPrintableBy<VarType, Printer>
{};

/// @brief Template alias for selecting a suitable printer class from the list of printers.
template<typename VarType, typename... Printers>
using SelectPrinter = ClassSelector<VarType, IsPrintableBy, Printers...>;





template<typename TupleType, typename Printer>
struct IsPrintableByWithTuple
{};


template<typename VarType, typename TraitsTuple, typename Printer>
struct IsPrintableByWithTuple<std::tuple<VarType, TraitsTuple>, Printer>
    : Impl::IsPrintableByWithTuple<VarType, TraitsTuple, Printer>
{};



template<typename VarType, typename TraitsTuple, typename Printer, typename... Printers>
using SelectPrinterTraitsTuple = ClassSelector<std::tuple<VarType, TraitsTuple>, IsPrintableByWithTuple, Printer, Printers...>;




template<typename T>
void VariableExport::exportVariable(std::ostream &ost, const T &var)
{
    SelectPrinter< T, PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                   PrintExportable, PrintTuple, PrintAny >::SelectedClass::print(ost, var);
}


template<typename T, typename ... TraitsTypes>
void VariableExport::exportVariable(std::ostream &ost, const T &var, const std::tuple<TraitsTypes...>& traitsTuple)
{
    SelectPrinterTraitsTuple< T, std::tuple<TraitsTypes...>,
                              PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                              PrintExportable, PrintTuple, PrintAny >::SelectedClass::print(ost, var, traitsTuple);
}


template<typename T>
void VariableExport::exportVariable(const T &var)
{
    SelectPrinter< T, PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                   PrintExportable, PrintTuple, PrintAny >::SelectedClass::print(var);
}


template<typename T, typename ... TraitsTypes>
void VariableExport::exportVariable(const T &var, const std::tuple<TraitsTypes...>& traitsTuple)
{
    SelectPrinterTraitsTuple< T, std::tuple<TraitsTypes...>,
                              PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                              PrintExportable, PrintTuple, PrintAny >::SelectedClass::print(var, traitsTuple);
}


#endif // VARIABLEEXPORT_H


