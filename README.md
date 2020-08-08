[![pipeline status](https://mmg-gitlab.fjfi.cvut.cz/gitlab/jakubec/UnstructuredMesh/badges/master/pipeline.svg)](https://mmg-gitlab.fjfi.cvut.cz/gitlab/jakubec/UnstructuredMesh/commits/master)

# GTMesh

The GTMesh is a C++ library utilizing modern C++ paradigms as template metaprogramming
and type traits. The aim of GTMesh is to provide an implementation working with an unstructured mesh
of any dimension and topology. Furthermore, the library provides additional tools developed
during the development of UnstructuredMesh and its functionalities. The minimal
required C++ standard is C++14.

The tools developed as part of GTMesh:
- [unstructured mesh](src/GTMesh/UnstructuredMesh/) with simple and user friendly inferface
- user friendly, simple and generic [debugging tool](src/GTMesh/Debug/)
- [Traits](src/GTMesh/Traits/), a tool describing C++ data structures and

## Basic Tutorial

Inicialization of the `UnstructuredMesh` object:

The `UnstructuredMesh` class has four template parameters, where the last is
variadic.

|Name|Type|Description|
|---|----|---|
|`Dimension`|`unsigned int`|Specifies the dimension of the stored mesh.|
|`IndexType`|`typename`|The type of the indices used in inner data structures.|
|`Real`|`typename`|The type of the coordinates of vertices.|
|`...Reserve`|`unsigned int`|Variadic parameter prescribing number of sub-elements of elements with dimension between 1 and Dimension-1. By default 0 is assumed, which implies dynamically allocated memory.|

The following code listing presents declaration of the `UnstructuredMesh`
object:
```c++
// 2D unstructured mesh size_t references and double precision coordinates
UnstructuredMesh<2, size_t, double> mesh2D;


// 3D unstructured mesh size_t references and double precision coordinates
// with maximum nuber of face sub-elements prescibed to 6
UnstructuredMesh<3, size_t, double, 6> mesh3D;

// with dynamically allocated sub-elements of faces,
// i.e. face could have any number of edges
UnstructuredMesh<3, size_t, double> mesh3D;


// 5D unstructured mesh size_t references and double precision coordinates
// maximum number of cell boundary sub-elements prescibed to 4
// maximum number of sub-elements of elements of dimension 3 prescibed to 4
// maximum number of sub-elements of elements of dimension 2 prescibed to 6
UnstructuredMesh<5, size_t, double, 4, 4, 6> mesh5D;
```
## Loading of mesh from a file
The file formats supproted by GTMesh are VTK (2D and 3D) and FPMA (3D only).
So far, only ASCII format of VTK file is supported.
```c++
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>

int main () {
    // create an instance of a 3D mesh with embedded memory 
    UnstructuredMesh<3, size_t, double, 6> mesh;
    
    // load mesh from vtk file
    auto meshReaderVTK = mesh.load("meshFile.vtk");
    
    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");
    return 0;
}
```

The reades objects contains additional information about the mesh such as
types of the cells. This information can be further used in the export of the
mesh.

## Mesh export
The following example loads a 3D mesh from FPMA file and exports it into VTK file.
```c++
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>

int main () {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    
    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");
    
    auto writer = mesh.write("meshFile1.vtk", *meshReaderFPMA, "My first exported mesh using GTMesh");
    
    return 0;
}
```

Note: Because the FPMA mesh format is able to store 3D mesh of arbitrary topology, the
mesh is exported tessellated into tetrahedrons. The object writer then contains
the mapping of original elements to the tessellated ones. This is used when mapping data
to the tessellated mesh.

## Mapping data to the mesh
For the purpose of mapping data to a mesh, GTMesh provides the
`MeshDataContainer` class. This data structure is able to allocate data separately
of the mesh. However, thanks to its interface it seems to the programmer 
as the data are allocated at the mesh itself. `MeshDataContainer` is able to 
store multiple data types to multiple dimensions of the mesh.
See the following example:
```c++
#include <GTMesh/Debug/Debug.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>

int main () {
    UnstructuredMesh<3, size_t, double, 6> mesh;

    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");

    // The flags container has internally 3 separate arrays (std::vector).
    // Each array is mapped to elements of different dimension (3, 2 and 1).
    // The data can be allocated to the mesh directly in the contructor,
    // also it is possible to initialize the allocated containers at once.
    MeshDataContainer<size_t, 3, 2, 1> flags(mesh, 1);

    MeshDataContainer<double, 3> cellsMeasures;
    // Memory allocation accoring to the mesh
    cellsMeasures.allocateData(mesh, 0);

    // Maps type int to the 3rd dimension and char to the 1st dimension
    MeshDataContainer<std::tuple<int, char>, 3,1> example(mesh, 0, 42);
    DBGVAR(example.getDataByPos<0>(), example.getDataByPos<1>());
    return 0;
}
```

There are two ways of addressing the inner arrays. The first is by its position
and the second is by dimension, the array is mapped to:
- `flags.getDataByPos<0>()` returns the array mapped as first (i.e., mapped to the 3rd dimension),
- `flags.getDataByDim<1>()` returns the array mapped to the 1st dimension (i.e., the array mapped as the 3rd).

Then, the data mapped to the first cell can be accessed by `flags.getDataByPos<0>()[0]`
or equivalently `flags.getDataByDim<3>()[0]`.
Finally, the data in `MeshDataContainer` can be adressed using the elements of
the mesh.
```c++
    for (auto& cell : mesh.getCells()){
        flags[cell] = cell.getIndex(); // Set the cell flag to the cell index
    }
```


## Mesh algorithms
The `UnstructuredMesh` class provides the following algorithms:
### Functions working with the mesh topology:
|Function name|Functionality|
|----|----|
|`apply`|Applies a given function to all connected elements of dimension one to elements of dimension two, where the dimensions are set as template parameters.|
|`color`|Generates a proper coloring of the mesh elements of dimension one connected by adjacency with elements of dimension two.|
|`connections`|Generates a MeshDataContainer containing the indexes of elements connected to another elements.|
|`neighbors`|Determines the neighborhood of elements based on 3 dimensions (start, connecting, target).|

### Functions calculating properties of mesh elements:
|Function name|Functionality|
|---|---|
|`computeElementCenters`|Calculates the centers of the mesh elements.|
|`initializeCenters`|Calculates the centers of the mesh elements and sets the centers of cells and faces up.|
|`computeElementMeasures`|Calculates measures of all elements except vertices in 2D and 3D meshes.|
|`computeFaceNormals`|Calculates normal vectors of faces in 2D and 3D meshes.|

Note: All these algorithms works in general dimension, e.g., 4D. In case of 3D
and 2D mesh, it is necessary for the centers of the mesh elements to be initialized
in order to the functions `computeElementMeasures` and `computeFaceNormals` work
properly.

The usage of the mesh algorithms is the following:
```c++
#include <GTMesh/Debug/Debug.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>

int main () {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    
    // Load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");
    
    mesh.apply<3, 0>(
        [&](size_t cellIndex, size_t vertIndex){
            DBGVAR(mesh.getVertices()[vertIndex]);
        }
    );
    
    // Determinte the vertices connected to cells
    auto connections3_0 = mesh.connections<3, 0>(); // MeshDataContainer<std::vector<size_t>, 3>
    
    // Vertices connected to the fisrst cell
    DBGVAR(connections3_0.getDataByPos<0>()[0]);
    
    // Determine a proper coloring of vertices with respect to connection by faces
    auto coloring0_2 = mesh.coloring<0, 2>();

    // Calculation of the mesh properties
    mesh.initializeCenters();
    auto measures = mesh.computeElementMeasures();
    
    auto normals = mesh.computeFaceNormals();
    return 0;
}
```

## Data export and import
The automatic data export and import is connected to the concept of [class traits](src/GTMesh/Traits/).
The only currently supported mesh format for data export is VTK. So far, only
data mapped to the cells can be exported.

```c++
#include <GTMesh/Debug/Debug.h>
#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataReader.h>
#include <GTMesh/UnstructuredMesh/MeshDataContainer/MeshDataIO/VTKMeshDataWriter.h>
#include <GTMesh/Traits/Traits.h>
#include <fstream>

struct cellData{
    double invVol;
    size_t index;
};

MAKE_ATTRIBUTE_TRAIT(cellData, invVol, index);

int main () {
    UnstructuredMesh<3, size_t, double, 6> mesh;

    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");

    MeshDataContainer<cellData, 3> meshData(mesh);

    mesh.initializeCenters();
    auto measures = mesh.computeElementMeasures();

    for (const auto& cell : mesh.getCells()){
        meshData[cell].invVol = 1.0 / measures[cell];
        meshData[cell].index = cell.getIndex();
    }

    std::ofstream oFile("meshFile-out.vtk");
    VTKMeshWriter<3, size_t, double> meshWriter;
    meshWriter.writeHeader(oFile, "mesh with data");
    meshWriter.writeToStream(oFile, mesh, meshReaderFPMA->getCellTypes());

    // append the data
    VTKMeshDataWriter<3>::writeToStream(oFile, meshData, meshWriter);

    // load the data back from meshFile-out
    UnstructuredMesh<3, size_t, double, 3> meshNew;
    auto meshReaderVTK = meshNew.load("meshFile-out.vtk");

    MeshDataContainer<cellData, 3> meshDataNew(meshNew);

    std::ifstream iFile("meshFile-out.vtk", std::ios::binary);

    VTKMeshDataReader<3, size_t>::readFromStream(iFile, meshDataNew);
    // Now the meshDataNew container contains the values stored in the meshFile-out.vtk
    DBGVAR(meshDataNew.getDataByDim<3>());
    return 0;
}
```