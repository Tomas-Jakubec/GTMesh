#ifndef DEBUG_H
#define DEBUG_H





#ifndef UNDEBUG
#include "../Macros/MacroForEach.h"
#include <iostream>
#include "HTMLLogger.h"
#include "CSVLogger.h"
#include "ConsoleLogger.h"
#include <stdexcept>
/*
** Macros intended for sending
** messages during program run
*/
namespace dbg {
    struct DBGStatics {
        static HtmlLogger HDBGLog;
        static CSVLogger CSVDBGLog;
    };

    HtmlLogger DBGStatics::HDBGLog("DBG.html");
    CSVLogger DBGStatics::CSVDBGLog("DBG.csv");
}

#define STRVAR(var) #var, var


#define DBGVAR(...) ConsoleLogger::writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND(condition, ...) if(condition) DBGVAR(__VA_ARGS__)

#define DBGMSG(...) ConsoleLogger::writeMessage("++", __LINE__, __FILE__, __VA_ARGS__)

#define DBGTRY(code)                        \
try{code;}                                   \
catch(const std::exception& e){              \
ConsoleLogger::writeMessage("!!", __LINE__, __FILE__, std::string("something went wrong in try block: ") + e.what());   \
abort();}

// Macros using html debug output
#define DBGVAR_HTML(...) dbg::DBGStatics::HDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND_HTML(condition, ...) if(condition) DBGVAR_HTML(__VA_ARGS__)

// Macros using csv debug output
#define DBGVAR_CSV(...) dbg::DBGStatics::CSVDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND_CSV(condition, ...) if(condition) DBGVAR_HTML(__VA_ARGS__)


#define DBGCHECK ConsoleLogger::writeMessage("--", __LINE__, __FILE__, "check line")


#else

#define DBG(comment)

#define DBGVARCOND(condition, ...)

#define DBGVAR(...)

#define DBGMSG(...)

#define DBGCHECK

#define DBGTRY(code) code

#define HTMLDBGVAR(...)

#define HTMLDBGCOND(condition, ...)

#define CSVDBGVAR(...)

#define CSVDBGCOND(condition, ...)
#endif //UNDEBUG

#endif // DEBUG_H
