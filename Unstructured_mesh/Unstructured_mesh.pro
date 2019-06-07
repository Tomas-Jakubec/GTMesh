TEMPLATE = app
CONFIG += console c++11
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
    mesh_element.h \
    unstructed_mesh_define.h \
    unstructuredmesh.h \
    vertex.h
