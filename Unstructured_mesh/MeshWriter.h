#ifndef MESHWRITER_H
#define MESHWRITER_H
#include "MeshNativeType.h"
#include "Vertex.h"
#include "MeshElement.h"

template<unsigned int MeshDimension>
class MeshWriter{


protected:


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

    template<typename IndexType, typename Real>
    struct HashData{
        IndexType nElem;
        Real totCoord;
    };

    template<typename IndexType, typename Real>
    union HashUni {
        HashData<IndexType, Real> data;
        char bytes[sizeof (HashData<IndexType, Real>) + 1] = {};
    };
public:
    using type = MeshNativeType<MeshDimension>;

    /**
     * @brief computeHash<HR>
     * Method calculating a hash from the mesh information to detect changes in the mesh.
     * The hash is calculated from number of elements and sum of vertices coordinates.
     * @param mesh input mesh to be hashed
     * @return hash code of type size_t
     */
    template<typename IndexType, typename Real, unsigned int ...Reserve>
    static size_t computeHash(MeshElements<MeshDimension, IndexType, Real, Reserve...>& mesh){

        // A vector of ones that simplifies the sum of coordinates
        Vertex<MeshDimension, Real> ones;
        for(unsigned int dim = 0; dim < MeshDimension; dim++){
            ones[dim] = 1.0;
        }

        Real totalVert = 0.0;
        for(IndexType i = 0; i < mesh.getVertices().size(); i++) {
            totalVert += mesh.getVertices().at(i) * ones;
        }

        IndexType numberOfElements = sumOfMeshElements<MeshDimension>::sum(mesh);

        HashUni<IndexType, Real> uni;
        uni.data = {numberOfElements, totalVert};
        std::hash<std::string> hasher;
        return hasher(uni.bytes);
        //return hasher(std::to_string(totalVert)+";"+std::to_string(numberOfElements));
    }
};


#endif // MESHWRITER_H
