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
