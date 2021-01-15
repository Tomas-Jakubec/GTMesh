#ifndef CONSOLELOGGER_H
#define CONSOLELOGGER_H
#include "VariableExport.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif


/**
 * @brief The ConsoleLogger class
 */
class ConsoleLogger {
    bool mColoredOutput = true; //!< If true then the output the variable names and their values are highlighted.
    bool mOnlyFilename = false; //!< If true then the filename is exported without its path.
    bool mPrintFunction = true; //!< If true then the name of the function the export is called from is printed.
    bool mPrintAggregated = true; //!< Prints the information about the origin of the message only once.
    std::string mIdentation = std::string("||     "); //!< The identation of the single exports when PrintAggregated is true.

public:
    /// If true then the output the variable names and their values are highlighted.
    void setColoredOutput(bool coloredOutput) { mColoredOutput = coloredOutput; }

    bool getColoredOutput() const { return mColoredOutput; }

    /// If true then the filename is exported without its path.
    void setOnlyFliename(bool onlyFilename) { mOnlyFilename = onlyFilename; }

    bool getOnlyFilename() const { return mOnlyFilename; }

    /// If true then the name of the function the export is called from is printed.
    void setPrintFunction(bool printFunction) { mPrintFunction = printFunction; }

    bool getPrintFunction() const { return mPrintFunction; }

    /// Prints the information about the origin of the message only once.
    void setPrintAggregated(bool printAggregated) { mPrintAggregated = printAggregated; }

    bool getPrintAggregated() const { return mPrintAggregated; }

     /// The identation of the single exports when PrintAggregated is true.
    void setIdentation(const std::string& identation) { mIdentation = identation; }

    std::string getIdentation() const { return mIdentation; }

    /// Resets the configuration to default
    void reset() {
        mColoredOutput = true; //!< If true then the output the variable names and their values are highlighted.
        mOnlyFilename = false; //!< If true then the filename is exported without its path.
        mPrintFunction = true; //!< If true then the name of the function the export is called from is printed.
        mPrintAggregated = true; //!< Prints the information about the origin of the message only once.
        mIdentation = std::string("...     "); //!< The identation of the single exports when PrintAggregated is true.
    }
public:

    template<typename VAR, typename ... REST>
    void writeVariable(int line, const char* funcName, const char* cppFile, const char* name,const VAR& value,const REST& ... rest){

        std::string fileName = getFileName(cppFile);

        if (mPrintAggregated) {
            printHeader(line, funcName, fileName.c_str());
            std::cerr << std::endl;
        }
        writeVar(line, funcName, fileName.c_str(), name, value, rest...);

    }

    template<typename MSGTYPE, typename... MSGTYPES>
    void writeMessage(const char* prefix, int line, const char* funcName, const char* sourceFile, const MSGTYPE& message, const MSGTYPES&... rest) {
        std::string fileName = getFileName(sourceFile);
        writeMsg(prefix, line, funcName, fileName.c_str(), message, rest...);
    }
public:
    template<typename MSGTYPE, typename... MSGTYPES>
    void writeMsg(const char* prefix, int line, const char* funcName, const char* sourceFile, const MSGTYPE& message, const MSGTYPES&... rest) {
        writeMsg(prefix, line, funcName, sourceFile, message);
        writeMsg(prefix, line, funcName, sourceFile, rest...);
    }




    template<typename MSGTYPE>
    void writeMsg(const char* prefix, int line, const char* funcName, const char* sourceFile, const MSGTYPE& message) {

        printHeader(line, funcName, sourceFile, prefix);
        std::cerr << "==> ";
        insertMessageColor();
        VariableExport::exportVariable(std::cerr, message);
        resetColor();
        std::cerr << " <==" << std::endl;

    }

private:
    template<typename VAR, typename ... REST>
    void writeVar(int line, const char* funcName, const char* cppFile, const char* name,const VAR& value,const REST& ... rest){

        writeVar(line, funcName, cppFile, name, value);
        writeVar(line, funcName, cppFile,  rest...);
    }

    template<typename VAR>
    void writeVar(int line, const char* funcName, const char* cppFile, const char* name,const VAR& value){

        if (mPrintAggregated) {
            std::cerr << mIdentation;
        } else {
            printHeader(line, funcName, cppFile);
            std::cerr << " ";
        }

        std::cerr << "[[ ";
        insertNameColor();
        std::cerr << name;
        resetColor();
        std::cerr << " ]] ==> " << std::flush;
        insertValueColor();
        VariableExport::exportVariable(std::cerr, value);
        resetColor();
        std::cerr << "\n";

    }

private:
    void insertNameColor() {
        if (mColoredOutput){
#ifdef _WIN32
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
#else
            std::cerr << "\033[0;33m";
#endif
        }
    }

    void resetColor() {
        if (mColoredOutput){
#ifdef _WIN32
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#else
            std::cerr << "\033[0m";
#endif
        }
    }

    void insertValueColor() {
        if (mColoredOutput){
#ifdef _WIN32
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_INTENSITY);
#else
            std::cerr << "\033[0;32m";
#endif
        }
    }

    void insertMessageColor() {
        if (mColoredOutput){
#ifdef _WIN32
            HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            SetConsoleTextAttribute(hConsole, 0x05 | FOREGROUND_INTENSITY);
#else
            std::cerr << "\033[0;31m";
#endif
        }
    }

    void printHeader(int line, const char* funcName, const char* cppFile, const char* prefix = "== ") {
        std::cerr << prefix << cppFile << ":" << line;
        if (mPrintFunction) {
            std::cerr << " <" << funcName << "> ";
        }
    }

    std::string getFileName(const char* cppFile){
        std::string res(cppFile);
        if (mOnlyFilename) {
#ifdef _WIN32
            size_t pos = res.find_last_of('\\');
#else
            size_t pos = res.find_last_of('/');
#endif
            if(pos != res.npos) {
                return res.substr(pos + 1);
            }
        }
        return res;
    }

public:

    template<typename MSGTYPE, typename... MSGTYPES>
    static void writeMessageC(const char* prefix, int line, const char* funcName, const char* sourceFile, const MSGTYPE& message, const MSGTYPES&... rest) {
        writeMessageC(prefix, funcName, line, sourceFile, message);
        writeMessageC(prefix, funcName, line, sourceFile, rest...);
    }


    template<typename MSGTYPE>
    static void writeMessageC(const char* prefix, int line, const char* funcName, const char* sourceFile, const MSGTYPE& message) {
        printf("%s %s: %i <%s> ==> ", prefix, sourceFile, line, funcName);

        VariableExport::exportVariable(message);

        printf(" <==\n");
    }



    template<typename VAR_NAME, typename VAR, typename ... REST>
    static void writeVarC(int line, const char* funcName, const char* cppFile, const VAR_NAME& name,const VAR& value,const REST& ... rest){

        writeVarC(line, funcName, cppFile, name, value);
        writeVarC(line, funcName, cppFile, rest...);
    }

    template<typename VAR>
    static void writeVarC(int line, const char* funcName, const char* cppFile, const char* name,const VAR& value){
        printf("== %s: %i <%s> [[ %s ]] ==> ", cppFile, line, funcName, name);

        VariableExport::exportVariable(value);
        printf("\n");
    }

};
#endif // CONSOLELOGGER_H
