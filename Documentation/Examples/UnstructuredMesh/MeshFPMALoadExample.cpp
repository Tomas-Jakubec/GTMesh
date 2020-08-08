#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>

int main () {
    UnstructuredMesh<3, size_t, double, 6> mesh;
    
    // load mesh from fpma file
    auto meshReaderFPMA = mesh.load("meshFile.fpma");
    
    auto writer = mesh.write("meshFile1.vtk", *meshReaderFPMA, "My first exported mesh using GTMesh");
    
    return 0;
}
