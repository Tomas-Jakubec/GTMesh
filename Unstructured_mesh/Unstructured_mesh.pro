TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    ../debug/debug.cpp

HEADERS += \
    cellboundaryconnection.h \
    cellconnection.h \
    ../debug/debug.h \
    ../debug/htmllogger.h \
    computationaly_significant_element.h \
    inline_array_operations.h \
    mesh_element.h \
    mesh_functions.h \
    unstructed_mesh_define.h \
    unstructuredmesh.h \
    vector.h \
    vertex.h
