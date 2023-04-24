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
