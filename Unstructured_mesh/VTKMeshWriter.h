#ifndef VTKMESHWRITER_H
#define VTKMESHWRITER_H

#include "MeshWriter.h"
#include "MeshElement.h"
#include "MeshDataContainer.h"
#include "MeshFunctions.h"
#include <map>
#include <ostream>
template<unsigned int MeshDimension, typename IndexType = size_t, typename Real = double>
class VTKMeshWriter : public MeshWriter<MeshDimension>{
public:
    VTKMeshWriter() = default;
    template<unsigned int ...Reserve>
    VTKMeshWriter(const MeshElements<MeshDimension, IndexType, Real, Reserve...>&){}
};



template<typename IndexType, typename Real>
class VTKMeshWriter<2, IndexType, Real> : public MeshWriter<2>{
    using writer = MeshWriter<2>;
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

    /**
     * @brief lastHash<HR>
     * The hash of the last written mesh.
     */
    size_t lastHash;
    /**
     * @brief cellVert<HR>
     * Vertices of all cells in correct order for vtk export.
     */
    MeshDataContainer<std::vector<IndexType>, 3> cellVert;

    /**
     * @brief totalNumberOfWrittenElements<HR>
     * Information required by VTK format.
     */
    IndexType totalNumberOfWrittenElements = 0;


    /**
     * @brief indexCell<HR>
     * This funcion stores indexes of vertices in correct order
     * in output vector verticesIndexed
     * @param mesh the structure of the mesh
     * @param face the face of the mesh to be indexed
     * @param cell the cell which are the vertices being indexed to
     * @param faceEdgeOri the orientation of the edges to the faces
     * @param verticesIndexed output vector of indexed vertices
     */
    template<unsigned int ...Reserve>
    void indexCell(MeshElements<2, IndexType, Real, Reserve...>& mesh,
                   typename MeshElements<2, IndexType, Real, Reserve...>::Cell& cell,
                   std::vector<IndexType>& verticesIndexed){


        IndexType tmpEdgeIndex = cell.getBoundaryElementIndex();

        IndexType nextVertex = INVALID_INDEX(IndexType);
        if (mesh.getEdges().at(tmpEdgeIndex).getCellLeftIndex() == cell.getIndex()){
            verticesIndexed.push_back(mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex());
            nextVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex();
        } else {
            verticesIndexed.push_back(mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex());
            nextVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex();
        }

        size_t lastWrittenEdge = tmpEdgeIndex;

        do {
            if ((lastWrittenEdge != tmpEdgeIndex)){
                if (mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex() == nextVertex) {

                    nextVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex();

                    verticesIndexed.push_back(mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex());

                    lastWrittenEdge = tmpEdgeIndex;

                } else if (mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex() == nextVertex){

                    nextVertex = mesh.getEdges().at(tmpEdgeIndex).getVertexBIndex();

                    verticesIndexed.push_back(mesh.getEdges().at(tmpEdgeIndex).getVertexAIndex());

                    lastWrittenEdge = tmpEdgeIndex;
                }


            }
            tmpEdgeIndex = mesh.getEdges().at(tmpEdgeIndex).getNextBElem(cell.getIndex());

        } while (nextVertex != verticesIndexed.at(0));

    }

    template<unsigned int ...Reserve>
    void indexMesh(MeshElements<2, IndexType, Real, Reserve...>& mesh){
        size_t curHash = writer::computeHash(mesh);

        // if the mesh is the same as it was, return
        if (lastHash == curHash){
            return;
        }
        lastHash = curHash;

        cellVert.template getDataByPos<0>().clear();

        std::vector<IndexType> vertIndex;
        for (auto& cell : mesh.getCells()) {
            vertIndex.clear();
            vertIndex.reserve(4);
            indexCell(mesh, cell, vertIndex);
            cellVert.template getDataByPos<0>().push_back(vertIndex);
        }
        cellVert.template getDataByPos<0>().shrink_to_fit();

        totalNumberOfWrittenElements = 0;
        for (auto verts : cellVert.template getDataByPos<0>()){
            totalNumberOfWrittenElements += verts.size() + 1;
        }


    }

    template<unsigned int ...Reserve>
    void writeToStream(std::ostream& ost,
                       MeshElements<2, IndexType, Real, Reserve...>& mesh,
                       MeshDataContainer<typename writer::type::ElementType, 2> cellTypes){

        indexMesh(mesh);
        // first write verices
        ost << "POINTS " << mesh.getVertices().size() <<
               " double" << std::endl;

        for(auto vert : mesh.getVertices()) {
            ost << vert[0] << " " << vert[1] << " 0.0\n";
        }
        ost << std::endl;


        ost << "CELLS " << mesh.getCells().size() << ' ' << totalNumberOfWrittenElements << std::endl;

        for (const auto& verts : cellVert.template getDataByPos<0>()){
            ost << verts.size() << ' ';
            for (IndexType i = 0; i < verts.size(); i++) {
                ost << verts[i];
                if (i < verts.size() - 1) {
                    ost << ' ';
                }
            }
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





template<typename IndexType, typename Real>
class VTKMeshWriter<3, IndexType, Real> : public MeshWriter<3>{
    using writer = MeshWriter<3>;
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

    // indexing of the mesh
    /**
     * @brief lastHash<HR>
     * The hash of the last written mesh.
     */
    size_t lastHash = 0;
    /**
     * @brief cellVert<HR>
     * Vertices of all cells in correct order for vtk export.
     */
    MeshDataContainer<std::vector<IndexType>, 3> cellVert;
    MeshDataContainer<typename writer::type::ElementType, 3> cellTypes;

    /**
     * @brief totalNumberOfWrittenElements<HR>
     * Information required by VTK format.
     */
    IndexType totalNumberOfWrittenElements = 0;


    std::vector<Vertex<3, Real>> appendedVertices;

    std::map<IndexType, IndexType> backwardCellIndexMapping;


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
                   typename MeshElements<3, IndexType, Real, Reserve...>::Cell& cell,
                   MeshDataContainer<std::vector<bool>, 2>& faceEdgeOri,
                   std::vector<IndexType>& verticesIndexed){


        IndexType startVertex = INVALID_INDEX(IndexType);
        IndexType nextVertex = INVALID_INDEX(IndexType);
        if (faceEdgeOri[face][0] == true && cell.getIndex() == face.getCellLeftIndex()){ // the edge is left to the face
            startVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexBIndex();
            nextVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexAIndex();
        } else {
            startVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexAIndex();
            nextVertex = mesh.getEdges().at(face.getSubelements()[0].index).getVertexBIndex();
        }

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

    }



private:
    /**
     * @brief indexPyramid<HR>
     * Stores indexes of vertices in output container vertWrit
     * @param mesh
     * @param cell the current proceeded cell
     * @param face the base face of the pyramid
     * @param faceEdgeOri orientation of edges to faces
     * @param vertWrit [out]
     */
    template<unsigned int ... Reserve>
    void indexPyramid(
            MeshElements<3, IndexType, Real, Reserve...>& mesh,
            typename MeshElements<3, IndexType, Real, Reserve...>::Cell& cell,
            typename MeshElements<3, IndexType, Real, Reserve...>::Face& face,
            MeshDataContainer<std::vector<bool>,2>& faceEdgeOri,
            std::vector<IndexType>& vertWrit
            ){
        indexFace(mesh,face, cell, faceEdgeOri, vertWrit);

        for (IndexType index : cellVert[cell]) {
            bool vertOK = true;
            for (IndexType i : vertWrit){
               if (i == index){
                   vertOK = false;
               }
            }
            if (vertOK){
                vertWrit.push_back(index);
            }
        }
    }




    /**
     * @brief indexLinearObject<HR>
     * This function stores indexes of vertices of cells which are made by cartesian product of
     * a line and the base face in output container vertWrit.
     * @param mesh
     * @param cell the current proceeded cell
     * @param face the base face of the pyramid
     * @param faceEdgeOri orientation of edges to faces
     * @param vertWrit [out]
     */
    template<unsigned int ... Reserve>
    void indexLinearObject(
            MeshElements<3, IndexType, Real, Reserve...>& mesh,
            typename MeshElements<3, IndexType, Real, Reserve...>::Cell& cell,
            typename MeshElements<3, IndexType, Real, Reserve...>::Face& face,
            MeshDataContainer<std::vector<bool>,2>& faceEdgeOri,
            std::vector<IndexType>& vertWrit
            ){
        indexFace(mesh, face, cell, faceEdgeOri, vertWrit);
        // write vertices of the oposite triangular side

        IndexType numVertBaseFace = vertWrit.size();

        IndexType index = 0;

        auto lambdaProc = [&mesh,&index,&vertWrit](IndexType, IndexType edgeIndex){
            auto& edge = mesh.getEdges().at(edgeIndex);

            if (edge.getVertexAIndex() == index){
                bool edgeOK = true;
                for (IndexType i : vertWrit){
                   if (edge.getVertexBIndex() == i){
                       edgeOK = false;
                   }
                }
                if(edgeOK){
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
                    vertWrit.push_back(edge.getVertexAIndex());
                }
            }
        };

        for (IndexType i = 0; i < numVertBaseFace; i++) {
            index = vertWrit[i];

            MeshRun<3,3,1,3,false, true>::run(
                mesh,
                cell.getIndex(),
                cell.getIndex(),
                lambdaProc
            );
        }
    }



public:

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
    void indexMesh(MeshElements<3, IndexType, Real, Reserve...>& mesh,
                   MeshDataContainer<typename writer::type::ElementType, 3> cellTypes){

        appendedVertices.clear();
        cellVert.template getDataByPos<0>().clear();
        backwardCellIndexMapping.clear();
        this->cellTypes.template getDataByPos<0>().clear();
DBGMSG("indexing mesh");
        // write cells of the mesh
        // prepare connections
        auto cellVert = MeshConnections<3,0>::connections(mesh);

        // prepare orientation for correct export
        // this is very expensive procedure
        auto faceEdgeOri = edgesOrientation(mesh);

        std::map<std::string, IndexType> appendedVertPos;



        std::vector<IndexType> vertWrit;
        for (auto cell : mesh.getCells()){
            vertWrit.clear();
            vertWrit.reserve(cellVert[cell].size());


            // switch the writing procedure according to cell type
            switch (cellTypes.template getDataByPos<0>().at(cell.getIndex())) {


            // Cell type TETRA
            case writer::type::ElementType::TETRA :{
                // every face is base face for TETRAHEDRON
                auto& face = mesh.getFaces().at(cell.getBoundaryElementIndex());

                indexPyramid(mesh, cell, face, faceEdgeOri, vertWrit);
                this->cellVert.template getDataByPos<0>().push_back(vertWrit);
                this->cellTypes.template getDataByPos<0>().push_back(cellTypes.at(cell));
            }break;



            // Cell type PYRAMID
            case writer::type::ElementType::PYRAMID :{

                // search for the base face
                typename MeshElements<3, IndexType, Real, Reserve...>::Face* face = nullptr;
                for (IndexType faceIndex : mesh.template getElement<3>(cell.getIndex()).getSubelements()){
                    if (mesh.template getElements<2>().at(faceIndex).getSubelements().getNumberOfSubElements() > 3){
                        face = &mesh.getFaces().at(faceIndex);
                    }
                }
                // index the pyramid object
                indexPyramid(mesh, cell, *face, faceEdgeOri, vertWrit);

                this->cellVert.template getDataByPos<0>().push_back(vertWrit);
                this->cellTypes.template getDataByPos<0>().push_back(cellTypes.at(cell));
            }break;



            // Cell type WEDGE
            case writer::type::ElementType::WEDGE :{
                // write vertices of one face
                typename MeshElements<3, IndexType, Real, Reserve...>::Face* face = nullptr;
                // search for the base face

                for (IndexType faceIndex : mesh.template getElement<3>(cell.getIndex()).getSubelements()){
                    if (mesh.template getElements<2>().at(faceIndex).getSubelements().getNumberOfSubElements() == 3){
                        face = &mesh.getFaces().at(faceIndex);
                        break;
                    }
                }

                indexLinearObject(mesh, cell, *face, faceEdgeOri, vertWrit);
                this->cellVert.template getDataByPos<0>().push_back(vertWrit);
                this->cellTypes.template getDataByPos<0>().push_back(cellTypes.at(cell));
            }break;



            // Cell type HEXAHEDRON
            case writer::type::ElementType::HEXAHEDRON :{
                // write vertices of one face
                auto& face = mesh.getFaces().at(cell.getBoundaryElementIndex());

                indexLinearObject(mesh, cell, face, faceEdgeOri, vertWrit);
                this->cellVert.template getDataByPos<0>().push_back(vertWrit);
                this->cellTypes.template getDataByPos<0>().push_back(cellTypes.at(cell));
            }break;


            // Defalt splitting of the volume
            default: {
                //throw std::runtime_error("it is not possible yet to write generic object into VTK");

                IndexType tmpFace = cell.getBoundaryElementIndex();
                IndexType cellCenterIndex = INVALID_INDEX(IndexType);
                std::string cellCenterKey = std::to_string(cell.getCenter()[0]) +
                        ";" + std::to_string(cell.getCenter()[1]) +
                        ";" + std::to_string(cell.getCenter()[2]);
                auto it = appendedVertPos.find(cellCenterKey);
                if (it == appendedVertPos.end()) {
                    appendedVertPos[cellCenterKey] = appendedVertices.size() + mesh.getVertices().size();
                    cellCenterIndex = appendedVertices.size() + mesh.getVertices().size();
                    appendedVertices.push_back(cell.getCenter());
                } else {
                    cellCenterIndex = it->second;
                }

                do {
                    auto& face = mesh.getFaces().at(tmpFace);

                    for (auto& sube : face.getSubelements()){
                        auto& edge = mesh.getEdges().at(sube.index);
                        vertWrit.clear();
                        vertWrit.reserve(4);
                        vertWrit.push_back(edge.getVertexAIndex());
                        vertWrit.push_back(edge.getVertexBIndex());
                        std::string faceCenterKey = std::to_string(face.getCenter()[0]) +
                                ";" + std::to_string(face.getCenter()[1]) +
                                ";" + std::to_string(face.getCenter()[2]);
                        IndexType faceCenterIndex = INVALID_INDEX(IndexType);
                        auto it = appendedVertPos.find(faceCenterKey);
                        if (it == appendedVertPos.end()) {
                            faceCenterIndex = appendedVertices.size() + mesh.getVertices().size();
                            appendedVertPos[faceCenterKey] = faceCenterIndex;
                            appendedVertices.push_back(face.getCenter());
                        } else {
                            faceCenterIndex = it->second;
                        }
                        vertWrit.push_back(faceCenterIndex);
                        vertWrit.push_back(cellCenterIndex);
                        backwardCellIndexMapping[this->cellVert.template getDataByPos<0>().size()] = cell.getIndex();
                        this->cellVert.template getDataByPos<0>().push_back(vertWrit);
                        this->cellTypes.template getDataByPos<0>().push_back(writer::type::ElementType::TETRA);
                    }

                    tmpFace = face.getNextBElem(cell.getIndex());
                } while (tmpFace != cell.getBoundaryElementIndex());

            }
            }
        }

        int cnt = 0;
        for (auto verts : this->cellVert.template getDataByPos<0>()){
            cnt += verts.size() + 1;
        }

        totalNumberOfWrittenElements = cnt;

    }


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
                       MeshElements<3, IndexType, Real, Reserve...>& mesh,
                       MeshDataContainer<typename writer::type::ElementType, 3> cellTypes){
        // create index of mesh if the mesh has changed
        size_t curHash = writer::computeHash(mesh);

        // if the mesh is not the same as it was,
        // then update the index
        if (lastHash != curHash){
            indexMesh(mesh, cellTypes);
        }
        lastHash = curHash;

        // first write verices
        ost << "POINTS " << mesh.getVertices().size() + appendedVertices.size() <<
               " double" << std::endl;

        for(auto vert : mesh.getVertices()) {
            ost << vert[0] << ' ' << vert[1] << ' ' << vert[2] <<"\n";
        }

        for(auto vert : appendedVertices) {
            ost << vert[0] << ' ' << vert[1] << ' ' << vert[2] <<"\n";
        }
        ost << std::endl;

        // write cells of the mesh
        ost << "CELLS " << cellVert.template getDataByPos<0>().size() << ' ' << totalNumberOfWrittenElements << std::endl;


        for (const auto& verts : cellVert.template getDataByPos<0>()){
            ost << verts.size() << ' ';
            for (IndexType i = 0; i < verts.size(); i++) {
                ost << verts[i];
                if (i < verts.size() - 1) {
                    ost << ' ';
                }
            }
            ost << std::endl;
        }

        ost << std::endl;

        ost << "CELL_TYPES " << this->cellTypes.template getDataByPos<0>().size() << std::endl;
        for (typename writer::type::ElementType type : this->cellTypes.template getDataByPos<0>()) {
            ost << TypeConversionTable.at(type) << "\n";
        }
        ost << std::endl;
    }
};


#endif // VTKMESHWRITER_H
