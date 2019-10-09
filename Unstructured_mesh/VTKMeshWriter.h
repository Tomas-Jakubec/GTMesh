#ifndef VTKMESHWRITER_H
#define VTKMESHWRITER_H

#include "MeshWriter.h"
#include "MeshElement.h"
#include "MeshDataContainer.h"
#include "MeshFunctions.h"
#include <map>
#include <ostream>
template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
class VTKMeshWriter : public MeshWriter<MeshDimension, IndexType, Real>{
public:
    VTKMeshWriter() = default;
    VTKMeshWriter(const MeshElements<MeshDimension, IndexType, Real, Reserve...>&){}
};



template<typename IndexType, typename Real, unsigned int... Reserve>
class VTKMeshWriter<2, IndexType, Real, Reserve...> : public MeshWriter<2, IndexType, Real>{
    using writer = MeshWriter<2, IndexType, Real>;
    std::map<typename writer::type::ElementType, int> TypeConversionTable{
        {writer::type::ElementType::LINE, 3},
        {writer::type::ElementType::TRIANGLE, 5},
        {writer::type::ElementType::QUAD, 9},
        {writer::type::ElementType::POLYGON, 7},
    };
public:
    void writeHeader(std::ostream& ost, const std::string& dataName) {
        ost << "# vtk DataFile Version 2.0\n" <<
               dataName << "\nASCII\nDATASET UNSTRUCTURED_GRID" << std::endl;
    }

    void writeToStream(std::ostream& ost,
                       MeshElements<2, IndexType, Real, Reserve...>& mesh,
                       MeshDataContainer<typename writer::type::ElementType, 2> cellTypes){
        // first write verices
        ost << "POINTS " << mesh.getVertices().size() <<
               " double" << std::endl;

        for(auto vert : mesh.getVertices()) {
            ost << vert[0] << " " << vert[1] << " 0.0\n";
        }
        ost << std::endl;

        auto cellVert = ::MeshConnections<2,0>::connections(mesh);
        int cnt = 0;
        for (auto verts : cellVert.template getDataByPos<0>()){
            cnt += verts.size() + 1;
        }
        ost << "CELLS " << mesh.getCells().size() << ' ' << cnt << std::endl;

        for (auto cell : mesh.getCells()){
            ost << cellVert.at(cell).size() << " ";

            size_t tmpEdgeIndex = cell.getBoundaryElementIndex();

            size_t lastWrittenVertex = INVALID_INDEX(IndexType);

            if (mesh.getEdges().at(tmpEdgeIndex).getCellLeftIndex() == cell.getIndex()){
                ost << mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex() << ' ';
                lastWrittenVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex();
            } else {
                ost << mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex() << ' ';
                lastWrittenVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex();
            }

            size_t lastWrittenEdge = tmpEdgeIndex;
            size_t verticesWritten = 1;
            do {
                if ((lastWrittenEdge != tmpEdgeIndex || lastWrittenEdge != INVALID_INDEX(IndexType) ) &&
                    (lastWrittenVertex == INVALID_INDEX(IndexType) || (lastWrittenVertex == mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex() || lastWrittenVertex == mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex()))) {
                    if (mesh.getEdges().at(tmpEdgeIndex).getCellLeftIndex() == cell.getIndex()) {

                        lastWrittenVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex();

                        ost << mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex();

                    } else {

                        lastWrittenVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex();

                        ost << mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex();

                    }

                    lastWrittenEdge = tmpEdgeIndex;

                    verticesWritten++;

                    ost << ' ';

                }
                tmpEdgeIndex = mesh.getEdges().at(tmpEdgeIndex).getNextBElem(cell.getIndex());

            } while (verticesWritten < cellVert.at(cell).size());
            ost << std::endl;
        }
        ost << std::endl;

        ost << "CELL_TYPES " << mesh.getCells().size() << std::endl;
        for (typename writer::type::ElementType type : cellTypes.template getDataByPos<0>()) {
            ost << TypeConversionTable.at(type) << "\n";
        }
        ost << std::endl;
    }
};





template<typename IndexType, typename Real, unsigned int... Reserve>
class VTKMeshWriter<3, IndexType, Real, Reserve...> : public MeshWriter<3, IndexType, Real>{
    using writer = MeshWriter<3, IndexType, Real>;
    std::map<typename writer::type::ElementType, int> TypeConversionTable{
        {writer::type::ElementType::TETRA, 10},
        {writer::type::ElementType::HEXAHEDRON, 12},
        {writer::type::ElementType::WEDGE, 13},
        {writer::type::ElementType::PYRAMID, 14},
    };

public:
    void writeHeader(std::ostream& ost, const std::string& dataName) {
        ost << "# vtk DataFile Version 2.0\n" <<
               dataName << "\nASCII\nDATASET UNSTRUCTURED_GRID" << std::endl;
    }


    std::vector<IndexType> writeFace(std::ostream& ost,
                   MeshElements<3, IndexType, Real, Reserve...>& mesh,
                   typename MeshElements<3, IndexType, Real, Reserve...>::Face& face,
                   MeshDataContainer<std::vector<bool>, 2>& faceEdgeOri){

        std::vector<IndexType> verticesWritten;
        verticesWritten.reserve(face.getSubelements().getNumberOfSubElements());

        IndexType startVertex = INVALID_INDEX(IndexType);
        IndexType nextVertex = INVALID_INDEX(IndexType);
        if (faceEdgeOri[face][0] == true){ // the edge is left to the face
            startVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexBIndex();
            ost << startVertex << ' ';
            nextVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexAIndex();
        } else {
            startVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexAIndex();
            ost << startVertex << ' ';
            nextVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexBIndex();
        }

        verticesWritten.push_back(startVertex);

        IndexType lastWrittenEdge = face.getSubelements()[0].index;
        while (startVertex != nextVertex){
            for (auto& sube : face.getSubelements()) {
                auto &edge = mesh.getEdges().at(sube.index);
                if (edge.getIndex() != lastWrittenEdge) {
                    if (edge.getVertexAIndex() == nextVertex) {
                        lastWrittenEdge = edge.getIndex();
                        ost << edge.getVertexAIndex() << ' ';
                        verticesWritten.push_back(edge.getVertexAIndex());
                        nextVertex = edge.getVertexBIndex();
                    } else if (edge.getVertexAIndex() == nextVertex) {
                        lastWrittenEdge = edge.getIndex();
                        ost << edge.getVertexBIndex() << ' ';
                        verticesWritten.push_back(edge.getVertexBIndex());
                        nextVertex = edge.getVertexAIndex();
                    }
                }
            }
        }
        return verticesWritten;
    }

    void writeToStream(std::ostream& ost,
                       MeshElements<3, IndexType, Real, Reserve...>& mesh,
                       MeshDataContainer<typename writer::type::ElementType, 3>& cellTypes){
        // first write verices
        ost << "POINTS " << mesh.getVertices().size() <<
               " double" << std::endl;

        for(auto vert : mesh.getVertices()) {
            ost << vert[0] << ' ' << vert[1] << ' ' << vert[2] <<"\n";
        }
        ost << std::endl;

        auto cellVert = MeshConnections<3,0>::connections(mesh);
        auto cellFace = MeshConnections<3,2>::connections(mesh);
        auto faceVert = MeshConnections<2,0>::connections(mesh);

        int cnt = 0;
        for (auto verts : cellVert.template getDataByPos<0>()){
            cnt += verts.size() + 1;
        }
        ost << "CELLS " << mesh.getCells().size() << ' ' << cnt << std::endl;

        auto faceEdgeOri = edgesOrientation(mesh);


        for (auto cell : mesh.getCells()){
            ost << cellVert.at(cell).size() << " ";



            switch (cellTypes.template getDataByPos<0>().at(cell.getIndex())) {

            case writer::type::ElementType::TETRA :{
                // write vertices of one face
                auto& face = mesh.getFaces().at(cell.getBoundaryElementIndex());

                std::vector<IndexType> vertWrit = writeFace(ost,mesh,face, faceEdgeOri);

                for (IndexType index : faceVert[face]) {
                    bool vertOK = true;
                    for (IndexType i : vertWrit){
                       if (i == index){
                           vertOK = false;
                       }
                    }
                    if (vertOK){
                        ost << index;
                        vertWrit.push_back(index);
                    }
                }
            }break;

            case writer::type::ElementType::PYRAMID :{
                // write vertices of one face
                typename MeshElements<3, IndexType, Real, Reserve...>::Face* face = nullptr;
                // search for the base face
                for (IndexType faceIndex : cellFace[cell]){
                    if (faceVert.template getDataByPos<0>().at(faceIndex).size() == 4){
                        face = &mesh.getFaces().at(faceIndex);
                    }
                }

                std::vector<IndexType> vertWrit = writeFace(ost, mesh, *face, faceEdgeOri);
                // write the last vertex
                for (IndexType index : faceVert.at(*face)) {
                    bool vertOK = true;
                    for (IndexType i : vertWrit){
                       if (i == index){
                           vertOK = false;
                       }
                    }
                    if (vertOK){
                        ost << index;
                        vertWrit.push_back(index);
                    }
                }
            }break;

            case writer::type::ElementType::WEDGE :{
                // write vertices of one face
                typename MeshElements<3, IndexType, Real, Reserve...>::Face* face = nullptr;
                // search for the base face
                for (IndexType faceIndex : cellFace[cell]){
                    if (faceVert.template getDataByPos<0>().at(faceIndex).size() == 3){
                        face = &mesh.getFaces().at(faceIndex);
                    }
                }

                std::vector<IndexType> vertWrit = writeFace(ost, mesh, *face, faceEdgeOri);
                // write vertices of the oposite triangular side
                for (IndexType index : vertWrit) {
                    MeshRun<3,3,1,3,false, true>::run(mesh, cell.getIndex(),cell.getIndex(),
                        [&ost,&mesh,&index,&vertWrit](IndexType, IndexType edgeIndex){
                        auto& edge = mesh.getEdges().at(edgeIndex);

                        if (edge.getVertexAIndex() == index){
                            bool edgeOK = true;
                            for (IndexType i : vertWrit){
                               if (edge.getVertexBIndex() == i){
                                   edgeOK = false;
                               }
                            }
                            if(edgeOK){
                                ost << edge.getVertexBIndex() << ' ';
                                vertWrit.push_back(edge.getVertexBIndex());
                            }
                        }

                        if (edge.getVertexBIndex() == index){
                            bool edgeOK = true;
                            for (IndexType i : vertWrit){
                               if (edge.getVertexAIndex() == i){
                                   edgeOK = false;
                               }
                            }
                            if(edgeOK){
                                ost << edge.getVertexAIndex() << ' ';
                                vertWrit.push_back(edge.getVertexAIndex());
                            }
                        }
                    }
                    );
                }
            }break;

            case writer::type::ElementType::HEXAHEDRON :{
                // write vertices of one face
                auto& face = mesh.getFaces().at(cell.getBoundaryElementIndex());

                std::vector<IndexType> vertWrit = writeFace(ost, mesh, face, faceEdgeOri);
                // write vertices of the oposite triangular side
                for (IndexType index : vertWrit) {
                    MeshRun<3,3,1,3,false, true>::run(mesh, cell.getIndex(),cell.getIndex(),
                        [&ost,&mesh,&index,&vertWrit](IndexType, IndexType edgeIndex){
                        auto& edge = mesh.getEdges().at(edgeIndex);

                        if (edge.getVertexAIndex() == index){
                            bool edgeOK = true;
                            for (IndexType i : vertWrit){
                               if (edge.getVertexBIndex() == i){
                                   edgeOK = false;
                               }
                            }
                            if(edgeOK){
                                ost << edge.getVertexBIndex() << ' ';
                                vertWrit.push_back(edge.getVertexBIndex());
                            }
                        }

                        if (edge.getVertexBIndex() == index){
                            bool edgeOK = true;
                            for (IndexType i : vertWrit){
                               if (edge.getVertexAIndex() == i){
                                   edgeOK = false;
                               }
                            }
                            if(edgeOK){
                                ost << edge.getVertexAIndex() << ' ';
                                vertWrit.push_back(edge.getVertexAIndex());
                            }
                        }
                    }
                    );
                }
            }break;
            default: throw std::runtime_error("it is not possible yet to write any object into VTK");
            }

            ost << std::endl;
        }
        ost << std::endl;

        ost << "CELL_TYPES" << mesh.getCells().size() << std::endl;
        for (typename writer::type::ElementType type : cellTypes.template getDataByPos<0>()) {
            ost << TypeConversionTable.at(type) << "\n";
        }
        ost << std::endl;
    }
};


#endif // VTKMESHWRITER_H
