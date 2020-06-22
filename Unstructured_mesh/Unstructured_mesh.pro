TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $$PWD/../src/

SOURCES += \
        main.cpp \

HEADERS += \
    ../src/GTMesh/Debug/JSONLogger.h \
    ../src/GTMesh/Macros/MacroForEach.h \
    ../src/GTMesh/Debug/CSVLogger.h \
    ../src/GTMesh/Debug/ConsoleLogger.h \
    ../src/GTMesh/Debug/Debug.h \
    ../src/GTMesh/Debug/HTMLLogger.h \
    ../src/GTMesh/Debug/VariableExport.h \
    ../src/GTMesh/NumericStaticArray/GramSchmidt.h \
    ../src/GTMesh/NumericStaticArray/InlineArrayOperations.h \
    ../src/GTMesh/Traits/CustomTypeTraits.h \
    ../src/GTMesh/Traits/MemberAccess/MemberAccess.h \
    ../src/GTMesh/Traits/TraitsAlgorithm/TraitsAlgorithm.h \
    ../src/GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataContainer.h \
    ../src/GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h \
    ../src/GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h \
    ../src/GTMesh/Singleton/Singleton.h \
    ../src/GTMesh/Traits/Traits.h \
    ../src/GTMesh/UnstructuredMesh/MeshElements/CellBoundaryConnection.h \
    ../src/GTMesh/UnstructuredMesh/MeshElements/CellConnection.h \
    ../src/GTMesh/UnstructuredMesh/MeshElements/ComputationalySignificantElement.h \
    ../src/GTMesh/UnstructuredMesh/MeshElements/MeshElement.h \
    ../src/GTMesh/UnstructuredMesh/MeshElements/MeshElements.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/ComputeCenters.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/ComputeMeasures.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/ComputeNormals.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/EdgesOrientation.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/MeshApply.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/MeshColoring.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/MeshConnections.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/MeshFunctions.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/MeshFunctionsDefine.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/CellsDistance.h \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/MeshNeighborhood.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshNativeType.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshReader/FPMAMeshReader.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshReader/MeshReader.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshReader/VTKMeshReader.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshWriter/FPMAMeshWriter.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshWriter/MeshWriter.h \
    ../src/GTMesh/UnstructuredMesh/MeshIO/MeshWriter/VTKMeshWriter.h \
    ../src/GTMesh/UnstructuredMesh/UnstructuredMesh.h \
    ../src/GTMesh/NumericStaticArray/Vector.h \
    ../src/GTMesh/NumericStaticArray/Vertex.h \
    ../src/GTMesh/UnstructuredMesh/UnstructuredMeshDefine.h \
    ../src/UnitTests/UnstructuredMesh/MeshSetup.h

DISTFILES += \
    ../src/GTMesh/Traits/README.md \
    ../src/UnitTests/Debug/DBGVAR_JSONTest.cpp \
    ../src/UnitTests/Debug/VariableExportTest.cpp \
    ../src/UnitTests/Traits/ArithmeticTraitsTest.cpp \
    ../src/UnitTests/Traits/TraitsTest.cpp \
    ../src/UnitTests/UnstructuredMesh/MeshDataContainerTest.cpp \
    ../src/UnitTests/UnstructuredMesh/UnstructuredMeshTest.cpp \
    ../README.md \
    ../src/GTMesh/Debug/README.md \
    ../src/GTMesh/Traits/TraitsAlgorithm/README.md \
    ../src/GTMesh/UnstructuredMesh/MeshFunctions/README.md\
    ../src/UnitTests/CMakeLists.txt \
    ../src/UnitTests/Debug/CMakeLists.txt \
    ../src/UnitTests/README.md \
    ../src/UnitTests/Traits/CMakeLists.txt \
    ../src/UnitTests/UnstructuredMesh/CMakeLists.txt
