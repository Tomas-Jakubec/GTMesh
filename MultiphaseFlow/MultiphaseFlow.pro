TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

SOURCES += \
        main.cpp \
    multiphaseflow.cpp

HEADERS += \
    multiphaseflow.h
