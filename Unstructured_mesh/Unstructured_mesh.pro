TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \

HEADERS += \
    ../src/Macros/MacroForEach.h \
    ../src/Debug/CSVLogger.h \
    ../src/Debug/ConsoleLogger.h \
    ../src/Debug/Debug.h \
    ../src/Debug/HTMLLogger.h \
    ../src/Debug/VariableExport.h \
    ../src/InlineArrayOperations.h \
    ../src/Traits/MemberApproach/MemberApproach.h \
    ../src/UnstructuredMesh/MeshDataContainer/MeshDataContainer.h \
    ../src/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h \
    ../src/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h \
    ../src/Singleton/Singleton.h \
    ../src/Traits/Traits.h \
    ../src/UnstructuredMesh/MeshElements/CellBoundaryConnection.h \
    ../src/UnstructuredMesh/MeshElements/CellConnection.h \
    ../src/UnstructuredMesh/MeshElements/ComputationalySignificantElement.h \
    ../src/UnstructuredMesh/MeshElements/MeshElement.h \
    ../src/UnstructuredMesh/MeshFunctions/MeshFunctions.h \
    ../src/UnstructuredMesh/MeshIO/MeshNativeType.h \
    ../src/UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h \
    ../src/UnstructuredMesh/MeshIO/MeshReader/MeshReader.h \
    ../src/UnstructuredMesh/MeshIO/MeshReader/VTKMeshReader.h \
    ../src/UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h \
    ../src/UnstructuredMesh/MeshIO/MeshWriter/MeshWriter.h \
    ../src/UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h \
    ../src/UnstructuredMesh/UnstructedMeshDefine.h \
    ../src/UnstructuredMesh/UnstructuredMesh.h \
    ../src/NumericStaticArray/Vector.h \
    ../src/NumericStaticArray/Vertex.h
