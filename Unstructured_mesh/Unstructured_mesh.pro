TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    ../debug/debug.cpp

HEADERS += \
    ../Macros/MacroForEach.h \
    ../debug/ConsoleLogger.h \
    ../debug/Debug.h \
    ../debug/HTMLLogger.h \
    ../debug/VariableExport.h \
    InlineArrayOperations.h \
    UnstructuredMesh/MeshDataContainer/MemberApproach.h \
    UnstructuredMesh/MeshDataContainer/MeshDataContainer.h \
    UnstructuredMesh/MeshDataContainer/Traits.h \
    UnstructuredMesh/MeshElements/CellBoundaryConnection.h \
    UnstructuredMesh/MeshElements/CellConnection.h \
    UnstructuredMesh/MeshElements/ComputationalySignificantElement.h \
    UnstructuredMesh/MeshElements/MeshElement.h \
    UnstructuredMesh/MeshFunctions/MeshFunctions.h \
    UnstructuredMesh/MeshIO/MeshNativeType.h \
    UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h \
    UnstructuredMesh/MeshIO/MeshReader/MeshReader.h \
    UnstructuredMesh/MeshIO/MeshReader/VTKMeshReader.h \
    UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h \
    UnstructuredMesh/MeshIO/MeshWriter/MeshWriter.h \
    UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h \
    UnstructuredMesh/UnstructedMeshDefine.h \
    UnstructuredMesh/UnstructuredMesh.h \
    Vector.h \
    Vertex.h
