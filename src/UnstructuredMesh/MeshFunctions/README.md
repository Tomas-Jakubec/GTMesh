## Mesh functions

The directory MeshFunctions includes header files implementing several important algorithms to work with an unstructured mesh.
One can include all the functions by including the MeshFunctions.h.

### Functions working with the mesh topology:
|class|static member function|functionality|
|----|----|----|
|MeshApply|apply|applies a given function to all connected elements of dimension one to elements of dimension two, where the dimensions
are set as template parameters|
|MeshColoring|color|generates a proper coloring of the mesh elements of dimension one connected by adjacency with elements of dimension two|

### List of functions calculating properties of mesh elements:
|Function name|functionality|
|------------|-------------|
|computeCenters|computes centers of mesh elements|
|cellsDistance|calculates the distances between neighboring cells centers|
|computeMeasures|calculates measures of all elements except vertices in 2D and 3D meshes|
|computeNormals|calculates normal vectors of faces in 2D and 3D meshes|
|edgesOrientation|determines whether an edge is left or right oriented to an adjacent face in 3D mesh|



Example how to calculate centers of mesh elements.
```c++
UnstructuredMesh<3, size_t, double> mesh;
MeshDataContainer<Vertex<3,double>, 1,2,3> centers = computeCenters(mesh);
// Or the same result can be obtained this way
auto centers = mesh.computeCenters;
```
