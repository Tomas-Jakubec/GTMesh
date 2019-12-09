#ifndef EDGESORIENTATION_H
#define EDGESORIENTATION_H
#include "../MeshElements/MeshElements.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/Vector.h"
#include "../../Debug/Debug.h"
#include "MeshApply.h"
#include "ComputeNormals.h"
#include <valarray>
#include <set>
#include <map>

template<typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<2, IndexType, Real, Reserve...>& mesh,
                typename MeshElements<2, IndexType, Real, Reserve...>::template ElementType<2>& face,
                typename MeshElements<2, IndexType, Real, Reserve...>::Edge& edge
                ) {

    Vertex<2, Real> AminC = mesh.getVertices().at(edge.getVertexAIndex()) - face.getCenter();
    Vertex<2, Real> BminC = mesh.getVertices().at(edge.getVertexBIndex()) - face.getCenter();



    Real res = AminC[0]*BminC[1]-BminC[0]*AminC[1];
    return res > 0;
    /*throw std::runtime_error("can not determine orientation of edge " +
                             std::to_string(edge.getIndex()) + " wrt face: " + std::to_string(face.getIndex()));*/
}

template<typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<3, IndexType, Real, Reserve...>& mesh,
                typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<2>& face,
                typename MeshElements<3, IndexType, Real, Reserve...>::Edge& edge,
                Vector<3, Real> faceNormal
                ) {

    Vertex<3, Real> AminC = mesh.getVertices().at(edge.getVertexAIndex()) - face.getCenter();
    Vertex<3, Real> BminC = mesh.getVertices().at(edge.getVertexBIndex()) - face.getCenter();

    Real res = Real(0);

    for (IndexType i = 0; i < 3; i++){
        IndexType ipo = (i+1)%(3);
        IndexType ipt = (i+2)%(3);
        res += AminC[i]*BminC[ipo]*faceNormal[ipt]-AminC[ipt]*BminC[ipo]*faceNormal[i];

    }
    return res > 0;
}

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh, IndexType faceIndex, IndexType edgeIndex) {

    typename MeshElements<MeshDimension, IndexType, Real, Reserve...>::Edge& edge = mesh.getEdges().at(edgeIndex);
    typename MeshElements<MeshDimension, IndexType, Real, Reserve...>::template ElementType<2>& face = mesh.template getElements<2>().at(faceIndex);

    return edgeIsLeft(mesh, face, edge);

}

template<typename IndexType, typename Real, unsigned int ...Reserve>
bool edgeIsLeft(MeshElements<3, IndexType, Real, Reserve...>& mesh, IndexType faceIndex, IndexType edgeIndex) {

    typename MeshElements<3, IndexType, Real, Reserve...>::Edge& edge = mesh.getEdges().at(edgeIndex);
    typename MeshElements<3, IndexType, Real, Reserve...>::template ElementType<2>& face = mesh.template getElements<2>().at(faceIndex);

    auto normals = ComputeFaceNormals<ComputationMethod::DEFAULT>(mesh);

    return edgeIsLeft(mesh, face, edge, normals[face]);

}

template<typename IndexType, typename Real, unsigned int ...Reserve>
MeshDataContainer<std::vector<bool>, 2> edgesOrientation(MeshElements<3, IndexType, Real, Reserve...>& mesh) {

    MeshDataContainer<std::vector<bool>, 2> orientations(mesh);
    auto normals = ComputeFaceNormals<ComputationMethod::DEFAULT>(mesh);

    for (auto& face : mesh.getFaces()) {
        orientations[face].resize(face.getSubelements().getNumberOfSubElements());
        for (IndexType i = 0; i < face.getSubelements().getNumberOfSubElements(); i++){
            typename MeshElements<3, IndexType, Real, Reserve...>::Edge& edge = mesh.getEdges().at(face.getSubelements()[i].index);

            orientations[face][i] = edgeIsLeft(mesh, face, edge, normals[face]);
        }
    }
    return orientations;
}


#endif // EDGESORIENTATION_H
