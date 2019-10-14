#ifndef MESHWRITER_H
#define MESHWRITER_H
#include "MeshNativeType.h"
#include "Vertex.h"
#include "MeshElement.h"

template<unsigned int MeshDimension>
class MeshWriter{


protected:
    /**
     * @brief The MeshHash struct<HR>
     * A container to store cumulative data about a mesh
     * to recognize changes in the mesh.
     * It uses number of elements and sum of vertices coordinations.
     */
    template<typename IndexType, typename Real>
    struct MeshHash {
        IndexType numberOfElements = 0;
        Real totalVert = 0;

        bool operator== (const MeshHash& rhs) const {
            return numberOfElements == rhs.numberOfElements &&
                    totalVert == rhs.totalVert;
        }

        bool operator!= (const MeshHash& rhs) const {
            return !((*this)==rhs);
        }
    };

private:
    template<unsigned int Dim, typename Dummy = void>
    struct sumOfMeshElements{
        template<typename IndexType, typename Real, unsigned int ...Reserve>
        static IndexType sum(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){
            return mesh.template getElements<Dim>().size() + sumOfMeshElements<Dim - 1>::sum(mesh);
        }
    };
    template<typename Dummy>
    struct sumOfMeshElements<0, Dummy>{
        template<typename IndexType, typename Real, unsigned int ...Reserve>
        static IndexType sum(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){
            return mesh.template getElements<0>().size();
        }
    };
public:
    using type = MeshNativeType<MeshDimension>;

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    static MeshHash<IndexType, Real> computeHash(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){
        MeshHash<IndexType, Real> res;
        // A vector of ones that simplifies the sum of coordinates
        Vertex<MeshDimension, Real> ones;
        for(unsigned int dim = 0; dim < MeshDimension; dim++){
            ones[dim] = 1.0;
        }

        for(IndexType i = 0; i < mesh.getVertices().size(); i++) {
            res.totalVert += mesh.getVertices().at(i) * ones;
        }

        res.numberOfElements = sumOfMeshElements<MeshDimension>::sum(mesh);
        return res;
    }
};


#endif // MESHWRITER_H
