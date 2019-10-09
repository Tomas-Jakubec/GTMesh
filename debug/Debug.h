#ifndef DEBUG_H
#define DEBUG_H



#define CONCATENATE(arg1, arg2)   CONCATENATE1(arg1, arg2)
#define CONCATENATE1(arg1, arg2)  CONCATENATE2(arg1, arg2)
#define CONCATENATE2(arg1, arg2)  arg1##arg2


#define FOR_EACH_1(what, x, ...) what(x)
#define FOR_EACH_2(what, x, ...)\
  what(x)\
  FOR_EACH_1(what,  __VA_ARGS__)
#define FOR_EACH_3(what, x, ...)\
  what(x)\
  FOR_EACH_2(what, __VA_ARGS__)
#define FOR_EACH_4(what, x, ...)\
  what(x)\
  FOR_EACH_3(what,  __VA_ARGS__)
#define FOR_EACH_5(what, x, ...)\
  what(x)\
 FOR_EACH_4(what,  __VA_ARGS__)
#define FOR_EACH_6(what, x, ...)\
  what(x)\
  FOR_EACH_5(what,  __VA_ARGS__)
#define FOR_EACH_7(what, x, ...)\
  what(x)\
  FOR_EACH_6(what,  __VA_ARGS__)
#define FOR_EACH_8(what, x, ...)\
  what(x)\
  FOR_EACH_7(what,  __VA_ARGS__)
#define FOR_EACH_9(what, x, ...)\
  what(x)\
  FOR_EACH_8(what,  __VA_ARGS__)
#define FOR_EACH_10(what, x, ...)\
  what(x)\
  FOR_EACH_9(what,  __VA_ARGS__)
#define FOR_EACH_11(what, x, ...)\
  what(x)\
  FOR_EACH_10(what,  __VA_ARGS__)
#define FOR_EACH_12(what, x, ...)\
  what(x)\
  FOR_EACH_11(what,  __VA_ARGS__)
#define FOR_EACH_13(what, x, ...)\
  what(x)\
  FOR_EACH_12(what,  __VA_ARGS__)
#define FOR_EACH_14(what, x, ...)\
  what(x)\
  FOR_EACH_13(what,  __VA_ARGS__)
#define FOR_EACH_15(what, x, ...)\
  what(x)\
  FOR_EACH_14(what,  __VA_ARGS__)
#define FOR_EACH_16(what, x, ...)\
  what(x)\
  FOR_EACH_15(what,  __VA_ARGS__)
#define FOR_EACH_17(what, x, ...)\
  what(x)\
  FOR_EACH_16(what,  __VA_ARGS__)



#define FOR_EACH_NARG(...) FOR_EACH_NARG_(__VA_ARGS__, FOR_EACH_RSEQ_N())
#define FOR_EACH_NARG_(...) FOR_EACH_ARG_N(__VA_ARGS__)
#define FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, N, ...) N
#define FOR_EACH_RSEQ_N() 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0

#define FOR_EACH_(N, what, ...) CONCATENATE(FOR_EACH_, N)(what, __VA_ARGS__)
#define FOR_EACH(what, ...) FOR_EACH_(FOR_EACH_NARG(__VA_ARGS__), what, __VA_ARGS__)


#define STRVAR(var) ,#var, var


#ifndef UNDEBUG
#include <iostream>
#include "HTMLLogger.h"
#include "ConsoleLogger.h"
#include <stdexcept>
/*
** Macros intended for sending
** messages during program run
*/

extern HtmlLogger HDBGLog;


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


#define DBGVAR(...) ConsoleLogger::writeVar(__LINE__, __FILE__ FOR_EACH(STRVAR, __VA_ARGS__))
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
#define HTMLDBGVAR(...) HDBGLog.writeVar(__LINE__, __FILE__ FOR_EACH(STRVAR, __VA_ARGS__))
#define HTMLDBGCOND(condition, ...) if(condition) HTMLDBGVAR(__VA_ARGS__)

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
