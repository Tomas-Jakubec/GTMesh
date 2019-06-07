#include <iostream>
#include "../debug/debug.h"
#include "unstructed_mesh_define.h"
#include "cellconnection.h"
#include "cellboundaryconnection.h"
#include "vertex.h"
#include "mesh_element.h"
#include "unstructuredmesh.h"
using namespace std;

int main()
{

    MeshElement<3,0,size_t, double> vrchol(1);
    vrchol.Vertex<3,double>::operator=( {123.0,156.0,684.0});


    MeshElement<3,0,size_t, double> vrchol2(2);
    vrchol2 = {123.0,156.0,684.0};

    auto a = vrchol2 += vrchol;

    vrchol2 *= 0.25;


    MeshElement<2,1,size_t, double> hrana(6);
    MeshElement<3,1,size_t, double> hrana3D(6);

    MeshElement<3,2,size_t,double> povrch(6);

    MeshElement<3,3,size_t, double> bunka(50);


    DBGVAR(vrchol, vrchol.NormEukleid(),vrchol2, vrchol2.GetIndex());

    UnstructuredMesh<3, size_t, double> mesh;

    auto& policko = mesh.GetElements<2>();

    policko.resize(100);

    policko.at(0).SetCellIndex(50);

    DBGVAR(mesh.GetElements<2>().at(0).GetCellLeftIndex());

    DBGVAR(INVALID_INDEX(size_t));
}
