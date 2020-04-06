## Mesh functions ##

The directory MeshFunctions includes header files implementing several important algorithms to work with an unstructured mesh.

|header name|functionality|
|------------|-------------|
|ComputeCenter.h|defines class computing centers of mesh elements|


Example how to calculate centers of mesh elements.
```c++
UnstructuredMesh<3, size_t, double> mesh;
MeshDataContainer<Vertex<3,double>, 1,2,3> centers = ComputeCenters(mesh);
```

</table>
