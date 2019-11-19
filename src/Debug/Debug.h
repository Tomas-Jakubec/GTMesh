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

#define DBG(comment)            \
std::cerr << "DBG in file "     \
<< __FILE__ << " at line " <<   \
__LINE__ << " says " << comment \
<< std::endl

#ifdef __linux__
#define SINGLE_DBGVAR(var)       \
std::cerr << "DBGVAR in file "   \
<< __FILE__ << " at line " <<    \
__LINE__ << " variable \033[0;33m" << #var \
<< "\033[0m has value of \033[0;31m" << var <<    \
"\033[0m" << std::endl
#elif _WIN32
#define SINGLE_DBGVAR(var)       \
std::cerr << "DBGVAR in file "   \
<< __FILE__ << " at line " <<    \
__LINE__ << " variable " << #var \
<< " has value of " << var << std::endl
#endif


#define SINGLE_DBGVAR_SC(var) SINGLE_DBGVAR(var);

#define DBGVAR(...) ConsoleLogger::writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND(condition, ...) if(condition) DBGVAR(__VA_ARGS__)

#ifdef __linux__
#define DBGMSG(message)       \
std::cerr << "DBGMSG contains \
message [[\033[0;32m" << message << "\033[0m]]"\
<< std::endl
#elif _WIN32
#define DBGMSG(message)       \
std::cerr << "DBGMSG contains \
message [[" << message << "]]"\
<< std::endl
#endif

#define DBGTRY(code)                        \
try{code;}                                   \
catch(const std::exception& e){              \
DBG("something went wrong in try block: " << e.what());   \
abort();}

// Macros using html debug output
#define HTMLDBGVAR(...) dbg::DBGStatics::HDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define HTMLDBGCOND(condition, ...) if(condition) HTMLDBGVAR(__VA_ARGS__)

// Macros using csv debug output
#define CSVDBGVAR(...) dbg::DBGStatics::CSVDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define CSVDBGCOND(condition, ...) if(condition) HTMLDBGVAR(__VA_ARGS__)


#define DBGCHECK DBG("check")

#else

#define SINGLE_DBGVAR_SC(var)

#define SINGLE_DBGVAR(var)

#define DBG(comment)

#define DBGVARCOND(condition, ...)

#define DBGVAR(...)

#define DBGMSG(msg)

#define DBGCHECK

#define DBGTRY(code) code

#define HTMLDBGVAR(...)

#define HTMLDBGCOND(condition, ...)
#endif //UNDEBUG

#endif // DEBUG_H
