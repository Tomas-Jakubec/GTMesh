#ifndef DEBUG_H
#define DEBUG_H

#ifndef UNDEBUG
#include "../Macros/MacroForEach.h"
#include <iostream>
#include "HTMLLogger.h"
#include "JSONLogger.h"
#include "ConsoleLogger.h"
#include "../Singleton/Singleton.h"
#include <stdexcept>
/*
** Macros intended for sending
** messages during program run
*/
namespace dbg {
    struct DBGStatics {
        HtmlLogger HDBGLog = HtmlLogger("DBG.html");

        JSONLogger JSONDBGLog = JSONLogger("DBG.json");

    };

    inline ConsoleLogger& getDebugConsoleLogger() {
        static ConsoleLogger consoleLogger = ConsoleLogger();
        return consoleLogger;
    }
}

#define STRVAR(var) #var, var


#define DBGVAR(...) dbg::getDebugConsoleLogger().writeVariable(__LINE__, __func__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND(condition, ...) if(condition) DBGVAR(__VA_ARGS__)

#define DBGVAR_STDIO(...) ConsoleLogger::writeVarC(__LINE__, __FILE__, __func__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND_STDIO(condition, ...) if(condition) DBGVAR_STDIO(__VA_ARGS__)

#define DBGMSG(...) dbg::getDebugConsoleLogger().writeMessage("++", __LINE__, __func__, __FILE__, __VA_ARGS__)

#define DBGTRY(code) \
try{code;} \
catch(const std::exception& e){ \
dbg::getDebugConsoleLogger().writeMessage("!!", __LINE__, __func__, __FILE__, std::string("something went wrong in try block: ") + e.what()); \
exit(1);}

// Macros using html debug output
#define DBGVAR_HTML(...) Singleton<dbg::DBGStatics>::getInstance().HDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND_HTML(condition, ...) if(condition) DBGVAR_HTML(__VA_ARGS__)

// Macros using csv debug output
#define DBGVAR_CSV(...) Singleton<dbg::DBGStatics>::getInstance().CSVDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND_CSV(condition, ...) if(condition) DBGVAR_CSV(__VA_ARGS__)

// Macros using json debug output
#define DBGVAR_JSON(...) Singleton<dbg::DBGStatics>::getInstance().JSONDBGLog.writeVar(__LINE__, __FILE__, FOR_EACH(STRVAR, __VA_ARGS__))
#define DBGVARCOND_JSON(condition, ...) if(condition) DBGVAR_JSON(__VA_ARGS__)



#define DBGCHECK dbg::getDebugConsoleLogger().writeMessage("--", __LINE__, __func__, __FILE__, "check line")
#define DBGCHECK_STDIO ConsoleLogger::writeMessageC("--", __LINE__, __func__, __FILE__, "check line")

#else

#define DBG(comment)

#define DBGVARCOND(condition, ...)

#define DBGVAR(...)

#define DBGMSG(...)

#define DBGCHECK

#define DBGTRY(code) code

#define DBGVAR_HTML(...)

#define DBGVARCOND_HTML(condition, ...)

#define DBGVAR_CSV(...)

#define DBGVARCOND_CSV(condition, ...)

#define DBGVAR_JSON(...)

#define DBGVARCOND_JSON(condition, ...)

#endif //UNDEBUG

#endif // DEBUG_H
