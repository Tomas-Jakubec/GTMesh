#ifndef GETCENTERS_H
#define GETCENTERS_H
#include "MeshFunctionsDefine.h"
#include "../MeshDataContainer/MeshDataContainer.h"
#include "../../NumericStaticArray/GramSchmidt.h"
#include <array>

namespace Impl {
template < unsigned int Dimension, unsigned int CurrentDimension, unsigned int MeshDimension >
struct GetCenters {
    template <typename IndexType, typename Real, unsigned int ...Reserve, unsigned int Dim = CurrentDimension, typename std::enable_if<Dim == 1, bool>::type = true>
    static void
    get( std::array< Vertex< MeshDimension, Real >, Dimension + 1 >& res,
         const MakeMeshDataContainer_t<Vertex< MeshDimension, Real >, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
         const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
         const IndexType& elementIndex ){


        (void) centers;
        res[Dimension - 1] = mesh.getVertices()[mesh.getEdges()[elementIndex].getVertexAIndex()];
        res[Dimension] = mesh.getVertices()[mesh.getEdges()[elementIndex].getVertexBIndex()];


    }

    template <typename IndexType, typename Real, unsigned int ...Reserve, unsigned int Dim = CurrentDimension, typename std::enable_if<Dim != MeshDimension && (Dim > 1), bool>::type = true>
    static void
    get( std::array< Vertex< MeshDimension, Real >, Dimension + 1 >& res,
         const MakeMeshDataContainer_t<Vertex< MeshDimension, Real >, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
         const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
         const IndexType& elementIndex ){

        res[Dimension - CurrentDimension] = centers.template getDataByDim<CurrentDimension>()[elementIndex];


        IndexType sube = mesh.template getElements<CurrentDimension>()[elementIndex].getSubelements()[0];
        GetCenters<Dimension, CurrentDimension - 1, MeshDimension>::get(res, centers, mesh, sube);

    }

    template <typename IndexType, typename Real, unsigned int ...Reserve, unsigned int Dim = CurrentDimension, typename std::enable_if<Dim == MeshDimension, bool>::type = true>
    static void
    get( std::array< Vertex< MeshDimension, Real >, Dimension + 1 >& res,
         const MakeMeshDataContainer_t<Vertex< MeshDimension, Real >, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
         const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
         const IndexType& elementIndex ){

        res[0] = centers.template getDataByDim<CurrentDimension>()[elementIndex];

        IndexType sube = mesh.template getElements<CurrentDimension>()[elementIndex].getBoundaryElementIndex();
        GetCenters<Dimension, CurrentDimension - 1, MeshDimension>::get(res, centers, mesh, sube);

        return res;

    }



};


/**
 * @brief This function returns an array of vertices
 * describing the "plane" the element is lying in.
 * The vertices are chosen as centers of subelements
 * (one per dimension) except edges. This choice prevents
 * obtaining linear dependent vertices. The processed element is
 * specified by its index @a elementIndex and dimension @a Dimension.
 * @param centers: centers of the mesh (can be obtained by the function calculateCenters)
 * @param mesh: mesh to be processed
 * @param elementIndex: index of the processed element
 */
template <unsigned int Dimension, unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
std::array< Vertex< MeshDimension, Real >, Dimension + 1 > getCenters( const MakeMeshDataContainer_t<Vertex< MeshDimension, Real >, make_custom_integer_sequence_t<unsigned int, 1, MeshDimension>>& centers,
                                                                       const MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh,
                                                                       IndexType elementIndex ) {
    std::array< Vertex< MeshDimension, Real >, Dimension + 1 > res;
    GetCenters<Dimension, Dimension, MeshDimension>::get(res, centers, mesh, elementIndex);
    return res;
}

}



#endif // GETCENTERS_H
