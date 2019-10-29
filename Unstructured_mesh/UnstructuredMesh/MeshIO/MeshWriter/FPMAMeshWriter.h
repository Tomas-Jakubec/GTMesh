#ifndef FPMAMESHWRITER_H
#define FPMAMESHWRITER_H

#include "MeshWriter.h"
#include "../../MeshElements/MeshElement.h"
#include "../../MeshDataContainer/MeshDataContainer.h"
#include "../../MeshFunctions/MeshFunctions.h"
#include <map>
#include <ostream>

template<unsigned int MeshDimension, typename IndexType = size_t, typename Real = double>
class FPMAMeshWriter : public MeshWriter<MeshDimension>{
public:
    FPMAMeshWriter() = default;
    template<unsigned int ...Reserve>
    FPMAMeshWriter(const MeshElements<MeshDimension, IndexType, Real, Reserve...>&){}
};


template<typename IndexType, typename Real>
class FPMAMeshWriter<3, IndexType, Real> : public MeshWriter<3>{

     using writer = MeshWriter<3>;
    size_t lastHash;
    /**
     * @brief faceVert<HR>
     * correctly erdered vertices for all faces
     */
    MeshDataContainer<std::vector<IndexType>, 2> faceVert;

    /**
     * @brief cellFace<HR>
     * original order of faces
     */
    MeshDataContainer<std::vector<IndexType>, 3> cellFace;

    /**
     * @brief indexFace<HR>
     * This funcion return indexes of vertices of a face in correct order
     * with respect to edge orientation
     * in output vector verticesIndexed
     * @param mesh the structure of the mesh
     * @param face the face of the mesh to be indexed
     * @param cell the cell which are the vertices being indexed to
     * @param faceEdgeOri the orientation of the edges to the faces
     * @param verticesIndexed output vector of indexed vertices
     */
    template<unsigned int ...Reserve>
    void indexFace(MeshElements<3, IndexType, Real, Reserve...>& mesh,
                   typename MeshElements<3, IndexType, Real, Reserve...>::Face& face,
                   std::vector<IndexType>& verticesIndexed){

        // export the face in "left" direction see VTK export
        IndexType startVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexBIndex();
        IndexType nextVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexAIndex();


        verticesIndexed.push_back(startVertex);

        IndexType lastWrittenEdge = face.getSubelements()[0].index;
        while (startVertex != nextVertex){
            for (auto& sube : face.getSubelements()) {
                auto &edge = mesh.getEdges().at(sube.index);

                if (edge.getIndex() != lastWrittenEdge) {
                    if (edge.getVertexAIndex() == nextVertex) {
                        lastWrittenEdge = edge.getIndex();
                        verticesIndexed.push_back(edge.getVertexAIndex());
                        nextVertex = edge.getVertexBIndex();
                    } else if (edge.getVertexBIndex() == nextVertex) {
                        lastWrittenEdge = edge.getIndex();
                        verticesIndexed.push_back(edge.getVertexBIndex());
                        nextVertex = edge.getVertexAIndex();
                    }
                }
            }
        }
        verticesIndexed.shrink_to_fit();
    }


    /**
     * @brief indexMesh<HR>
     * This function creates vector of indexes of vertices for each cell in order suitable
     * for VTK output. Moreover in case of elements type different from
     * any VTK cell type, the cell is split into tetrahedrons. The tetrahedrons are
     * made for each edge of every faces. The tetrahedrons is made of both edges vertices
     * and cell and face center.
     * @param mesh Mesh to be indexed
     * @param cellTypes Vector of known cell types. If a cell type is not known then a method
     * of splitting into tetrahedrons is used.
     */
    template<unsigned int ...Reserve>
    void indexMesh(MeshElements<3, IndexType, Real, Reserve...>& mesh){


        faceVert.template getDataByPos<0>().clear();
        faceVert.alocateData(mesh);
DBGMSG("indexing mesh");
        // write cells of the mesh
        // prepare connections
        auto cellVert = MeshConnections<3,0>::connections(mesh);

        for (typename MeshElements<3, IndexType, Real, Reserve...>::Face& face : mesh.getFaces()){
            indexFace(mesh, face, faceVert.at(face));
        }

        cellFace = MeshConnections<3,2, Order::ORDER_ORIGINAL>::connections(mesh);
    }

public:
    /**
     * @brief writeToStream
     * Exports the mesh to the output stream in VTK format.
     * If the mesh is written for the first time, it is indexed.
     * @see indexMesh
     * @param ost
     * @param mesh
     * @param cellTypes the types of cells in NativeType format
     */
    template<unsigned int ...Reserve>
    void writeToStream(std::ostream& ost,
                       MeshElements<3, IndexType, Real, Reserve...>& mesh){
        // create index of mesh if the mesh has changed
        size_t curHash = writer::computeHash(mesh);

        // if the mesh is not the same as it was,
        // then update the index
        if (lastHash != curHash){
            indexMesh(mesh);
        }
        lastHash = curHash;

        // first write verices
        ost << mesh.getVertices().size() << std::endl;

        for(auto& vert : mesh.getVertices()) {
            ost << vert[0] << ' ' << vert[1] << ' ' << vert[2] <<"\n";
        }

        ost << std::endl;

        // write faces
        ost << mesh.getFaces().size() << std::endl;
        for(auto& face : mesh.getFaces()) {
            ost << faceVert.at(face).size() << ' ';
            for (IndexType& vertIndex : faceVert.at(face)) {
                ost << vertIndex << ' ';
            }
            ost << std::endl;
        }

        // write cells of the mesh
        ost << mesh.getCells().size() << std::endl;
        for(auto& cell : mesh.getCells()) {
            ost << cellFace.at(cell).size() << ' ';
            for (IndexType& faceIndex : cellFace.at(cell)) {
                ost << faceIndex << ' ';
            }
            ost << std::endl;
        }
        ost << 0;
    }
};






#endif // FPMAMESHWRITER_H
