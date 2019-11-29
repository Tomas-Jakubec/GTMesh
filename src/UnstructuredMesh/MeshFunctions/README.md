## Mesh functions ##

The directory MeshFunctions includes headder files implementing several important algorithms to work with an unstructured mesh.
<table>
<tr>
    <th>headder name</th><th>functionality</th>
</tr>
<tr>
    <ti>ComputeCenter.h</ti><ti>defines class computing centers of mesh elements</ti>
</tr>

Example how to calculate centers of mesh elements.
```c++
UnstructuredMesh<3, size_t, double> mesh;
MeshDataContainer<Vertex<3,double>, 1,2,3> centers = ComputeCenters(mesh);
```

</table>
