#ifndef HTMLLOGGER_H
#define HTMLLOGGER_H
#include "VariableExport.h"
#include <iostream>
#include <fstream>
#include <string>

/**
 * @brief The HtmlLogger class
 * has opened file to write log
 * if it opens existing file overwrites it.
 * in destructor it closes the file.
 */
class HtmlLogger {

    std::string logFileName = "";
    std::ofstream* logFile;
    int groupIndex = 0;

    inline void writeScript();
    inline void writeEndScritpt();
public:
    HtmlLogger(){logFile = nullptr;}

    /**
     * @brief HtmlLogger
     * opens specified file for writing
     * @param fileName
     *
     */
    HtmlLogger(const char* fileName)
        :HtmlLogger(){
        logFileName = fileName;
    }

    ~HtmlLogger(){
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


        writeScript();
    }

    void destroy(){
        if(logFile){
            writeEndScritpt();
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

        (*logFile)  << "[" << groupIndex << ", " << line << ", '" << cppFile << "', '" << name << "', ";
        VariableExport<>::exportVariable(*logFile,value);
        (*logFile) << "],\n";
        writeVar(line, cppFile, rest...);
    }


    template<typename VAR_NAME, typename VAR>
    void writeVar(int line, std::string& cppFile, VAR_NAME name, VAR value){

        // create file if not logFile is nullptr
        if (!logFile) create(logFileName.c_str());

        (*logFile)  << "[" << groupIndex << ", " << line << ", '" << cppFile << "', '" << name << "', ";
        VariableExport<>::exportVariable(*logFile,value);
        (*logFile) << "],\n";
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




void HtmlLogger::writeScript()
{
    (*logFile) << R"(<html>
<head>
</head>
<body>
<input type="number" id ="selected_line" value="-1"/>
<table id="selector">
<tr><button onclick="addFilter() " >add filter</button></tr>
</table><br/>
<button onclick="copyData(); showData() ">show all data</button>
<button onclick="filterData(); showData() ">filter</button>
<hr/>
<div id = "listDiv"></div>
<script type="text/javascript">
var data;
var dataToShow = new Array();
var numberOfSelectors = 0;
var apliedRule;
function copyData() {
    dataToShow = data;
}

function addFilter() {
    var sel = document.getElementById("selector");

    sel.innerHTML += '<tr><td>select: <select id ="sO' + numberOfSelectors.toString(10) + '"><option>row</option><option>group</option></select></td><td>var name:<input id ="vN' + numberOfSelectors.toString(10) + '" type="text"/></td><td><select id ="cO' + numberOfSelectors.toString(10) + '"><option>less</option><option>greater</option><option>equal</option><option>between</option></select></td><td><input type="number" step="1e-6" id ="firstVal' + numberOfSelectors.toString(10) + '" value="0"/></td><td> and: <input type="number" step="1e-6" id ="secondVal' + numberOfSelectors.toString(10) + '" value="0"/></td></tr>';

    numberOfSelectors++;
}

function checkData(varData) {
    if (numberOfSelectors > 0){
        for (var i = 0; i < numberOfSelectors; i++) {
            var selectOption = document.getElementById("sO" + i.toString(10)).selectedIndex;
            var varName = document.getElementById("vN" + i.toString(10)).value;
            var option = document.getElementById("cO" + i.toString(10)).selectedIndex;
            var val1 = document.getElementById("firstVal" + i.toString(10)).value;
            var val2 = document.getElementById("secondVal" + i.toString(10)).value;

            if (varData[3] == varName || varData[3] == "") {
                switch (option){

                case 0: {if(varData[4] < val1) {apliedRule = selectOption; return true;}} break;
                case 1: {if(varData[4] > val1) {apliedRule = selectOption; return true;}} break;
                case 2: {if(varData[4] == val1) {apliedRule = selectOption; return true;}} break;
                case 3: {if(varData[4] > val1 && varData[4] < val2) {apliedRule = selectOption; return true;}} break;

                }
            }
        }
    } else {
        return true;
    }

    return false;
}

function filterData() {

    dataToShow = [];

    var restrLine = document.getElementById("selected_line").value;


    for(var i = 0; i < data.length; i++){
        if (data[i][1] == restrLine || restrLine < 0) {

            if (checkData(data[i])){
                switch(apliedRule){
                    case 0 : dataToShow.push(data[i]); break;
                    case 1 : {
                        var index = i;
                        var groupIndex = data[i][0];
                        for(var j = i; j >= 0 && data[j][0] == groupIndex; j--){
                            index = j;
                        }
                        while(index < data.length && data[index][0] == groupIndex) {
                            dataToShow.push(data[index]);
                            index++;
                        }
                        i = index - 1;
                    } break;
                }
            }
        }
    }
}


function showData(){
    lD = document.getElementById("listDiv");

    var tmp = "<table style=\"position: absolute; width:90%; left:5%\"><tr style=\"background-color:lightgrey\"><th>GroupIndex</th><th>Line</th><th>File</th><th>var name</th><th>value</th></tr>";
    var odd = false;
    var oddGroup = false;
    if (dataToShow.length > 0) {
        var curGroup = dataToShow[0][0];
    }
    for(var i = 0; i < dataToShow.length; i++){

        if (curGroup != dataToShow[i][0]) {
            curGroup = dataToShow[i][0];
            oddGroup = !oddGroup;
        }

        if (odd) {
            tmp += "<tr style=\"background-color:" + (oddGroup ? "#ccffcc" : "#ffeecc") + "\"><td>";
            odd = false;
        } else {
            tmp += "<tr style=\"background-color:" + (oddGroup ? "#66ff66" : "#ffcc66") + "\"><td>";
            odd = true;
        }

        tmp += dataToShow[i][0] + "</td><td>" + dataToShow[i][1] + "</td><td>" + dataToShow[i][2] + "</td><td>" + dataToShow[i][3] + "</td><td>" + JSON.stringify(dataToShow[i][4]) + "</td></tr>";

    }
    tmp += "</table>";
    lD.innerHTML = tmp;

}

</script>


<script type = "text/javascript">

data = [)";
}

void HtmlLogger::writeEndScritpt()
{
    std::cerr << "writing end of the script" << std::endl;
    *logFile << R"(];
copyData();
showData();
</script>
</body>
</html>
)";
}



#endif // HTMLLOGGER_H
