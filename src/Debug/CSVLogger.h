#ifndef CSVLOGGER_H
#define CSVLOGGER_H
#include <string>
#include <fstream>
#include "VariableExport.h"
/**
 * @brief The HtmlLogger class
 * has opened file to write log
 * if it opens existing file overwrites it.
 * in destructor it closes the file.
 */
class CSVLogger {

    std::string logFileName = "";
    std::ofstream* logFile;
    int groupIndex = 0;

    void writeHeadder() {
        (*logFile) << "Group Index;Line;File;Expression;Value\n";
    }
public:
    CSVLogger(){logFile = nullptr;}

    /**
     * @brief HtmlLogger
     * opens specified file for writing
     * @param fileName
     *
     */
    CSVLogger(const char* fileName)
        :CSVLogger(){
        logFileName = fileName;
    }

    ~CSVLogger(){
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
        logFile = new std::ofstream(fileName);


        writeHeadder();
    }

    void destroy(){
        if(logFile){
            logFile->close();
            delete logFile;
        }
    }



    template<typename VAR_NAME, typename VAR, typename ... REST>
    void writeVar(VAR_NAME name, VAR value, REST ... rest){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());

        (*logFile) << name << ", " << value;
        logFile->flush();
        writeVar(rest...);
    }


    template<typename VAR_NAME, typename VAR>
    void writeVar(VAR_NAME name, VAR value){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());

        (*logFile) << name << ", " << value;
        logFile->flush();
    }


    template<typename VAR_NAME, typename VAR, typename ... REST>
    void writeVar(int line, std::string& cppFile, VAR_NAME name, VAR value, REST ... rest){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());

        (*logFile)  <<  groupIndex << ';' << line << ";\"" << cppFile << "\";\"" << name << "\";";
        VariableExport::_writeWar(*logFile,value);
        (*logFile) << "\n";
        writeVar(line, cppFile, rest...);
    }


    template<typename VAR_NAME, typename VAR>
    void writeVar(int line, std::string& cppFile, VAR_NAME name, VAR value){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());

        (*logFile)  <<  groupIndex << ';' << line << ";\"" << cppFile << "\";\"" << name << "\";";
        VariableExport::_writeWar(*logFile,value);
        (*logFile) << "\n";

        groupIndex++;
    }


    template<typename VAR_NAME, typename VAR, typename ... REST>
    void writeVar(int line, const char* cppFile, VAR_NAME name, VAR value, REST ... rest){

        std::string file;
        int i = 0;
        while(cppFile[i] != '\0'){
            if (cppFile[i] == '\\'){
                file += "\\\\";
            } else {
                file += cppFile[i];
            }
            i++;
        }
        writeVar(line, file, name, value,  rest...);
    }


    template<typename VAR_NAME, typename VAR>
    void writeVar(int line, const char* cppFile, VAR_NAME name, VAR value){
        std::string file;
        int i = 0;
        while(cppFile[i] != '\0'){
            if (cppFile[i] == '\\'){
                file += "\\\\";
            } else {
                file += cppFile[i];
            }
            i++;
        }
        writeVar(line, file, name, value);
    }

};


#endif // CSVLOGGER_H
