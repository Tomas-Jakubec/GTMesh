#ifndef JSONLOGGER_H
#define JSONLOGGER_H
#include <string>
#include <fstream>
#include "VariableExport.h"
#include <iostream>
#include <iomanip>
/**
 * @brief The JSONLogger class
 * has opened file to write log
 * if it opens existing file overwrites it.
 * in destructor it closes the file.
 */
class JSONLogger {

    std::string logFileName = "";
    std::ofstream* logFile;
    int groupIndex = 0;
    bool firstWrite = true;

public:
    JSONLogger(){logFile = nullptr;}

    /**
     * @brief JSONLogger
     * opens specified file for writing
     * @param fileName
     *
     */
    JSONLogger(const char* fileName)
        : JSONLogger(){
        logFileName = fileName;
    }

    ~JSONLogger(){
        destroy();
    }

    /**
     * @brief Create
     * opens specified file and writes headder and script in html
     * @param fileName
     */
    void create(const char* fileName){
        if (logFile){
            logFile->close();
            delete logFile;
        }
        firstWrite = true;
        logFile = new std::ofstream(fileName);
        (*logFile) << "{\n\"logs\":[\n";
    }

    void destroy(){
        if(logFile){
            (*logFile) << "\n]\n}";
            logFile->flush();

            logFile->close();
            delete logFile;
        }
    }



    template<typename VAR, typename ... REST>
    void writeVar(const char* name, const VAR& value, REST ... rest){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());

        writeVar(name, value);
        writeVar(rest...);
    }


    template<typename VAR>
    void writeVar(const char* name, const VAR& value){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());
        if (!firstWrite) {
            (*logFile) << ",\n";
        } else {
            firstWrite = false;
        }

        (*logFile) << "\t{\n\t\t\"expr\" : \""<< std::quoted( name ) <<
                      "\",\n\t\t\"data\" : " ;
        VariableExport<>::exportVariable(*logFile, value);
        (*logFile) << "\n}";


        logFile->flush();
    }


    template<typename VAR, typename ... REST>
    void writeVar(int line, const std::string& cppFile, const char* name, const VAR& value, const REST& ... rest){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());
        if (!firstWrite) {
            (*logFile) << ",\n";
        } else {
            firstWrite = false;
        }
        (*logFile) << "\t{" <<
                      "\n\t\t\"gInd\" : " << groupIndex << "," <<
                      "\n\t\t\"file\" : \"" << cppFile << "\"," <<
                      "\n\t\t\"line\" : " << line << "," <<
                      "\n\t\t\"expr\" : \""<< std::quoted( name ) << "\"," <<
                      "\n\t\t\"data\" : ";
        VariableExport<>::exportVariable(*logFile, value);
        (*logFile) << "\n\t}";

        writeVar(line, cppFile, rest...);

    }


    template<typename VAR>
    void writeVar(int line, const std::string& cppFile, const char* name, const VAR& value){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());
        if (!firstWrite) {
            (*logFile) << ",\n";
        } else {
            firstWrite = false;
        }
        (*logFile) << "\t{" <<
                      "\n\t\t\"gInd\" : " << groupIndex << "," <<
                      "\n\t\t\"file\" : \"" << cppFile << "\"," <<
                      "\n\t\t\"line\" : " << line << "," <<
                      "\n\t\t\"expr\" : \""<< std::quoted( name ) << "\"," <<
                      "\n\t\t\"data\" : ";
        VariableExport<>::exportVariable(*logFile, value);

        (*logFile) << "\n\t}";

        groupIndex++;
    }


    template< typename VAR, typename ... REST >
    void writeVar(int line, const char* cppFile, const char* name, const VAR& value,const REST& ... rest){
        writeVar(line, doubleBackSlash(cppFile), name, value,  rest...);
    }



    std::string doubleBackSlash(const char * str){
        std::string res = str;
        size_t pos = 0;
        while((pos = res.find("\\", pos + 2)) != res.npos){
            res.replace(pos, 1, "\\\\");
        }
        return res;
    }
};

#endif // JSONLOGGER_H
