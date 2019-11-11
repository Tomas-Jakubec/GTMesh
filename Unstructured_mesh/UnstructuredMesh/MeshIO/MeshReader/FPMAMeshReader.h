#ifndef FPMAMESHREADER_H
#define FPMAMESHREADER_H
#include "MeshReader.h"
#include "../../MeshElements/MeshElement.h"
#include <iostream>
#include <unordered_map>

template <unsigned int MeshDimension>
class FPMAMeshReader : public MeshReader<MeshDimension>{
};

template <>
class FPMAMeshReader<3> : public MeshReader<3> {
public:




    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadVertices(std::istream& ist,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numVert;
        ist >> numVert;

        mesh.getVertices().resize(numVert);

        for (IndexType i = 0; i < numVert; i++) {

            typename MeshElements<3, IndexType, Real, Reserve...>::Vertex& vert = mesh.getVertices().at(i);
            vert.setIndex(i);
            ist >> vert[0];
            ist >> vert[1];
            ist >> vert[2];
        }
    }


    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadFaces(std::istream& ist,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        // map of constructed edges
        std::unordered_map<std::string, IndexType> edges;

        IndexType numFace;
        ist >> numFace;

        mesh.getFaces().resize(numFace);

        // read face verts
        for (IndexType faceIndex = 0; faceIndex < numFace; faceIndex++) {
            mesh.getFaces().at(faceIndex).setIndex(faceIndex);
            IndexType numVert;
            ist >> numVert;

            std::vector<IndexType> vertices(numVert);
            for(IndexType j = 0; j < numVert; j++){
                ist >> vertices.at(j);
            }


            for(IndexType j = 0; j < numVert; j++){

                IndexType iA = vertices.at(j), iB = vertices.at((j+1)%numVert);
                std::string edgeKey = iA < iB ? std::to_string(iA) +";"+ std::to_string(iB) : std::to_string(iB) +";"+ std::to_string(iA);
                typename std::unordered_map<std::string, IndexType>::iterator edgeIt = edges.find(edgeKey);

                IndexType edgeIndex = IndexType();

                if (edgeIt == edges.end()){

                    edgeIndex = mesh.getEdges().size();
                    mesh.getEdges().push_back({});
                    mesh.getEdges().at(edgeIndex).setVertexAIndex(iA);
                    mesh.getEdges().at(edgeIndex).setVertexBIndex(iB);
                    mesh.getEdges().at(edgeIndex).setIndex(edgeIndex);

                    edges[edgeKey] = edgeIndex;
                } else {
                    edgeIndex = edgeIt->second;
                }
                try {
                    mesh.getFaces().at(faceIndex).getSubelements().addSubelement(edgeIndex);
                } catch (std::runtime_error& err) {
                    throw std::runtime_error(std::string("The number of edges has overflew the prealocated memory ") +
                                             std::to_string(mesh.getFaces().at(faceIndex).getSubelements().size()) +
                                             " in face number " + std::to_string(faceIndex) + " while adding edge number " +
                                             std::to_string(edgeIndex) + " ( vertA : " +std::to_string(iA) +
                                             "vertB : " + std::to_string(iB) + ")\n" + err.what());
                }
            }

        }
        mesh.getEdges().shrink_to_fit();
    }


    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadCells(std::istream& ist,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;

        mesh.getCells().resize(numCells);
        // read cells
        for (IndexType cellIndex = 0; cellIndex < numCells; cellIndex++) {
            mesh.getCells().at(cellIndex).setIndex(cellIndex);

            // read cell faces
            IndexType numFaces;
            ist >> numFaces;


            IndexType prevFace = INVALID_INDEX(IndexType);
            for(IndexType j = 0; j < numFaces; j++){
                IndexType faceIndex;
                ist >> faceIndex;

                if (j == 0){
                    mesh.getCells().at(cellIndex).setBoundaryElementIndex(faceIndex);
                }

                if (prevFace != INVALID_INDEX(IndexType)){
                    mesh.getFaces().at(prevFace).setNextBElem(faceIndex, cellIndex);

                }

                if (j == numFaces -1) {
                    mesh.getFaces().at(faceIndex).setNextBElem(mesh.getCells().at(cellIndex).getBoundaryElementIndex(), cellIndex);
                }
                prevFace = faceIndex;
            }


        }

    }



    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadFromStream(std::istream& ist,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        ist.seekg(ist.beg);

        loadVertices(ist, mesh);

        loadFaces(ist, mesh);

        loadCells(ist, mesh);

    }


};


#endif // FPMAMESHREADER_H
