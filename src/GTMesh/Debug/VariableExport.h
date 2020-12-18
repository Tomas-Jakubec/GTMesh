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




enum VARIABLE_EXPORT_METHOD {
    VARIABLE_EXPORT_METHOD_OSTREAM,
    VARIABLE_EXPORT_METHOD_STDIO
};

template <VARIABLE_EXPORT_METHOD target = VARIABLE_EXPORT_METHOD::VARIABLE_EXPORT_METHOD_OSTREAM>
struct VariableExport {


    template <typename T>
    static void exportVariable(std::ostream& ost, const T& var);


    template <typename T, typename ... TraitsTypes>
    static void exportVariable(std::ostream& ost, const T& var, const std::tuple<TraitsTypes...>&);
};

template <>
struct VariableExport<VARIABLE_EXPORT_METHOD_STDIO>
{

    template <typename T>
    static void exportVariable(const T& var);


    template <typename T, typename ... TraitsTypes>
    static void exportVariable(const T& var, const std::tuple<TraitsTypes...>&);

};


/*
 * Printers
 */
#include "Printers/PrintAny.h"
#include "Printers/PrintCustom.h"
#include "Printers/PrintExportable.h"
#include "Printers/PrintText.h"
#include "Printers/PrintIterable.h"
#include "Printers/PrintIndexable.h"
#include "Printers/PrintTuple.h"
#include "Printers/PrintTraitedClass.h"

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

template<typename VarType, typename Printer, typename... Printers>
struct SelectPrinter : public std::conditional_t<IsPrintableBy<VarType, Printer>::value,
                                                 SelectPrinter<VarType, Printer>,
                                                 SelectPrinter<VarType, Printers...>>
{};

template<typename VarType, typename Printer>
struct SelectPrinter<VarType, Printer>
{
    using PrinterType = Printer;
};

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

template<typename VarType, typename TraitsTuple, typename Printer, typename... Printers>
struct SelectPrinterTraitsTuple
    : public std::conditional_t<IsPrintableByWithTuple<VarType, TraitsTuple, Printer>::value,
                                SelectPrinterTraitsTuple<VarType, TraitsTuple, Printer>,
                                SelectPrinterTraitsTuple<VarType, TraitsTuple, Printers...>>
{};

template<typename VarType, typename TraitsTuple, typename Printer>
struct SelectPrinterTraitsTuple<VarType, TraitsTuple, Printer>
{
    using PrinterType = Printer;
};

template<>
template<typename T>
void VariableExport<VARIABLE_EXPORT_METHOD_OSTREAM>::exportVariable(std::ostream &ost, const T &var)
{
    SelectPrinter< T, PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                   PrintExportable, PrintTuple, PrintAny >::PrinterType::print(ost, var);
}

template<>
template<typename T, typename ... TraitsTypes>
void VariableExport<VARIABLE_EXPORT_METHOD_OSTREAM>::exportVariable(std::ostream &ost, const T &var, const std::tuple<TraitsTypes...>& traitsTuple)
{
    SelectPrinterTraitsTuple< T, std::tuple<TraitsTypes...>,
                              PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                              PrintExportable, PrintTuple, PrintAny >::PrinterType::print(ost, var, traitsTuple);
}


template<typename T>
void VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(const T &var)
{
    SelectPrinter< T, PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                   PrintExportable, PrintTuple, PrintAny >::PrinterType::print(var);
}


template<typename T, typename ... TraitsTypes>
void VariableExport<VARIABLE_EXPORT_METHOD_STDIO>::exportVariable(const T &var, const std::tuple<TraitsTypes...>& traitsTuple)
{
    SelectPrinterTraitsTuple< T, std::tuple<TraitsTypes...>,
                              PrintCustom, PrintText, PrintIterable, PrintIndexable, PrintTraitedClass,
                              PrintExportable, PrintTuple, PrintAny >::PrinterType::print(var, traitsTuple);
}


#endif // VARIABLEEXPORT_H


