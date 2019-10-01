TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    ../debug/debug.cpp

HEADERS += \
    CellBoundaryConnection.h \
    CellConnection.h \
    ComputationalySignificantElement.h \
    InlineArrayOperations.h \
    MeshDataContainer.h \
    MeshElement.h \
    MeshFunctions.h \
    MeshNativeType.h \
    MeshReader.h \
    MeshWriter.h \
    UnstructedMeshDefine.h \
    UnstructuredMesh.h \
    VTKMeshReader.h \
    ../debug/debug.h \
    ../debug/htmllogger.h \
    VTKMeshWriter.h \
    Vector.h \
    Vertex.h
