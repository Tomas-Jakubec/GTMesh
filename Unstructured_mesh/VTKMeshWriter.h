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

        auto cellVert = temp1::MeshConnections<2,0>::connections(mesh);
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
class VTKMeshWriter<3, IndexType, Real, Reserve...> : public MeshWriter<2, IndexType, Real>{
    using writer = MeshWriter<2, IndexType, Real>;
    std::map<typename writer::type::ElementType, int> TypeConversionTable{
        {writer::type::ElementType::TETRA, 10},
        {writer::type::ElementType::HEXAGON, 12},
        {writer::type::ElementType::WEDGE, 13},
        {writer::type::ElementType::PYRAMID, 14},
    };

public:
    void writeHeader(std::ostream& ost, const std::string& dataName) {
        ost << "# vtk DataFile Version 2.0\n" <<
               dataName << "\nASCII\nDATASET UNSTRUCTURED_GRID" << std::endl;
    }

    void writeToStream(std::ostream& ost,
                       MeshElements<3, IndexType, Real, Reserve...>& mesh,
                       MeshDataContainer<typename writer::type::ElementType, 3> cellTypes){
        // first write verices
        ost << "POINTS " << mesh.getVertices().size() <<
               " double" << std::endl;

        for(auto vert : mesh.getVertices()) {
            ost << vert[0] << ' ' << vert[1] << ' ' << vert[2] <<"\n";
        }
        ost << std::endl;

        auto cellVert = temp1::MeshConnections<3,0>::connections(mesh);
        auto cellFace = temp1::MeshConnections<3,2>::connections(mesh);
        auto faceVert = temp1::MeshConnections<2,0>::connections(mesh);

        int cnt = 0;
        for (auto verts : cellVert.template getDataByPos<0>()){
            cnt += verts.size() + 1;
        }
        ost << "CELLS " << mesh.getCells().size() << ' ' << cnt << std::endl;

        for (auto cell : mesh.getCells()){
            ost << cellVert.at(cell).size() << " ";



            switch (cellTypes.template getDataByPos<0>().at(cell.getIndex())) {

            case writer::type::ElementType::TETRA :{
                // write vertices of one face
                auto& face = mesh.getFaces(*(cellFace.at(cell).begin()));

                int vertWritten = 0;
                IndexType nextVertex = INVALID_INDEX(IndexType);
                if (face.getSubelements[0].isLeft){
                    ost << mesh.getEdges().at(face.getSubelements[0].index).getVetexAIndex();
                    nextVertex = mesh.getEdges().at(face.getSubelements[0].index).getVertexBIndex() << ' ';
                } else {
                    ost << mesh.getEdges().at(face.getSubelements[0].index).getVetexBIndex();
                    nextVertex = mesh.getEdges().at(face.getSubelements[0].index).getVertexAIndex() << ' ';
                }

                IndexType lastWrittenEdge = face.getSubelements[0].index;
                while (vertWritten != faceVert.at(face).size()){
                    for (auto& sube : face.getSubelements()) {
                        auto &edge = mesh.getEdges().at(sube.index);
                        if (edge.getIndex() != lastWrittenEdge) {
                            if (edge.getVertexAIndex() == nextVertex) {
                                lastWrittenEdge = edge.getIndex();
                                ost << edge.getVertexAIndex() << ' ';
                                nextVertex = edge.getVertexAIndex();
                            }
                        }
                    }
                }
            }break;


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
