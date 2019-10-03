#ifndef VTKMESHREADER_H
#define VTKMESHREADER_H

#include "MeshReader.h"
#include "MeshDataContainer.h"
#include "MeshElement.h"
#include <istream>
#include <string>
#include <unordered_map>
#include <map>

template<unsigned int MeshDimension, typename IndexType, typename Real, unsigned int ...Reserve>
class VTKMeshReader : public MeshReader<MeshDimension, IndexType, Real>{
public:
    VTKMeshReader() = default;
    VTKMeshReader(const MeshElements<MeshDimension, IndexType, Real, Reserve...>&){}
};



template<typename IndexType, typename Real, unsigned int... Reserve>
class VTKMeshReader<2, IndexType, Real, Reserve...> : public MeshReader<2, IndexType, Real>{
    using reader = MeshReader<2, IndexType, Real>;
    std::map<int, typename reader::type::ElementType> TypeConversionTable{
        {3, reader::type::ElementType::LINE},
        {5, reader::type::ElementType::TRIANGLE},
        {8, reader::type::ElementType::QUAD},
        {9, reader::type::ElementType::QUAD},
        {7, reader::type::ElementType::POLYGON},
    };

    std::unordered_map<std::string, IndexType> edges;

    MeshDataContainer<typename reader::type::ElementType, 2> cellTypes;
    // file indexing
    //

    //
    //MeshDataContainer<IndexType>
public:
    VTKMeshReader() = default;
    VTKMeshReader(const MeshElements<2, IndexType, Real, Reserve...>&){}

    MeshDataContainer<typename reader::type::ElementType, 2> getCellTypes() {
        return cellTypes;
    }

    void loadPoints(std::istream& ist, MeshElements<2, IndexType, Real, Reserve...>& mesh){
        IndexType numPoints;
        ist >> numPoints;
        Real dummy = 0;
        mesh.getVertices().resize(numPoints);
        ist.ignore(20, '\n');
        for (IndexType vertIndex = 0; vertIndex < numPoints; vertIndex++) {
            mesh.getVertices().at(vertIndex).setIndex(vertIndex);
            ist >> mesh.getVertices().at(vertIndex)[0];
            ist >> mesh.getVertices().at(vertIndex)[1];
            ist >> dummy;
        }
    }

    void loadCells(std::istream& ist, MeshElements<2, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;
        mesh.getCells().resize(numCells);
        // Skip the total number of numbers
        ist.ignore(50, '\n');
        for (IndexType cellIndex = 0; cellIndex < numCells; cellIndex++) {
            mesh.getCells().at(cellIndex).setIndex(cellIndex);
            IndexType numVert;
            ist >> numVert;

            std::vector<IndexType> vertices(numVert);
            for(IndexType j = 0; j < numVert; j++){
                ist >> vertices.at(j);
            }

            IndexType prevEdge = INVALID_INDEX(IndexType);
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

                    mesh.getEdges().at(edgeIndex).setCellLeftIndex(cellIndex);
                } else {
                    edgeIndex = edgeIt->second;
                    mesh.getEdges().at(edgeIt->second).setCellRightIndex(cellIndex);
                }

                if (prevEdge != INVALID_INDEX(IndexType)){
                    mesh.getEdges().at(prevEdge).setNextBElem(edgeIndex, cellIndex);
                }

                if (j == 0){
                    mesh.getCells().at(cellIndex).setBoundaryElementIndex(edgeIndex);
                }
                if (j == numVert - 1) {
                    mesh.getEdges().at(edgeIndex).setNextBElem(mesh.getCells().at(cellIndex).getBoundaryElementIndex(), cellIndex);
                }
                prevEdge = edgeIndex;
            }

        }
    }


    void loadCellTypes(std::istream& ist, MeshElements<2, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;
        cellTypes.alocateData(mesh);
        for (IndexType i = 0; i < numCells; i++) {
            int vtkType = 0;
            ist >> vtkType;
            typename std::map<int, typename reader::type::ElementType>::iterator typeIt = TypeConversionTable.find(vtkType);
            if (typeIt != TypeConversionTable.end()){
                cellTypes.template getDataByPos<0>().at(i) = typeIt->second;
            } else {
                std::runtime_error("unsuported cell type");
            }
        }
    }


    void loadFromStream(std::istream& ist,MeshElements<2, IndexType, Real, Reserve...>& mesh){
        ist.seekg(ist.beg);
        // Ignore first row "# vtk DataFile Version 2.0"
        ist.ignore(1024, '\n');
        // Ignore name of the data set
        ist.ignore(1024, '\n');
        // ASCII or BINARY
        std::string buf;
        std::getline(ist, buf);
        if (buf != "ASCII"){
            throw std::runtime_error("ASCII expected but got " + buf);
        }

        ist >> buf;
        if (buf != "DATASET"){
            throw std::runtime_error("the keyword DATASET expected");
        }

        ist >> buf;
        if (buf != "UNSTRUCTURED_GRID"){
            throw std::runtime_error("only unstructured grid is supported but got " + buf);
        }


        ist >> buf;
        if (buf == "POINTS") {
            loadPoints(ist, mesh);
        }

        ist >> buf;
        if (buf == "CELLS") {
            loadCells(ist, mesh);
        }

        ist >> buf;
        if (buf == "CELL_TYPES") {
            loadCellTypes(ist, mesh);
        }

    }

    MeshElements<2, IndexType, Real, Reserve...> loadFromStream(std::istream& ist){
        MeshElements<2, IndexType, Real, Reserve...> resultMesh;
        loadFromStream(ist, resultMesh);
        return resultMesh;
    }
};


template<typename IndexType, typename Real, unsigned int... Reserve>
class VTKMeshReader<3, IndexType, Real, Reserve...> : public MeshReader<3, IndexType, Real>{
    using reader = MeshReader<3, IndexType, Real>;
    std::map<int, typename reader::type::ElementType> TypeConversionTable{
        {10, reader::type::ElementType::TETRA},
        {11, reader::type::ElementType::HEXAHEDRON},
        {12, reader::type::ElementType::HEXAHEDRON},
        {13, reader::type::ElementType::WEDGE},
        {14, reader::type::ElementType::PYRAMID},
    };

    std::unordered_map<std::string, IndexType> edges;
    std::unordered_map<std::string, IndexType> faces;

    MeshDataContainer<typename reader::type::ElementType, 3> cellTypes;
    // file indexing
    //

    //
    //MeshDataContainer<IndexType>
public:
    VTKMeshReader() = default;
    VTKMeshReader(const MeshElements<3, IndexType, Real, Reserve...>&){}

    MeshDataContainer<typename reader::type::ElementType, 3> getCellTypes() {
        return cellTypes;
    }

    void loadPoints(std::istream& ist, MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numPoints;
        ist >> numPoints;
        Real dummy = 0;
        mesh.getVertices().resize(numPoints);
        ist.ignore(20, '\n');
        for (IndexType vertIndex = 0; vertIndex < numPoints; vertIndex++) {
            mesh.getVertices().at(vertIndex).setIndex(vertIndex);
            ist >> mesh.getVertices().at(vertIndex)[0];
            ist >> mesh.getVertices().at(vertIndex)[1];
            ist >> mesh.getVertices().at(vertIndex)[2];
        }
    }

    void loadCells(std::istream& ist, MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;
        mesh.getCells().resize(numCells);
        // Skip the total number of numbers
        ist.ignore(50, '\n');
        for (IndexType cellIndex = 0; cellIndex < numCells; cellIndex++) {
            mesh.getCells().at(cellIndex).setIndex(cellIndex);
            IndexType numVert;
            ist >> numVert;

            std::vector<IndexType> vertices(numVert);
            for(IndexType j = 0; j < numVert; j++){
                ist >> vertices.at(j);
            }

            IndexType prevEdge = INVALID_INDEX(IndexType);
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

                    mesh.getEdges().at(edgeIndex).setCellLeftIndex(cellIndex);
                } else {
                    edgeIndex = edgeIt->second;
                    mesh.getEdges().at(edgeIt->second).setCellRightIndex(cellIndex);
                }

                if (prevEdge != INVALID_INDEX(IndexType)){
                    mesh.getEdges().at(prevEdge).setNextBElem(edgeIndex, cellIndex);
                }

                if (j == 0){
                    mesh.getCells().at(cellIndex).setBoundaryElementIndex(edgeIndex);
                }
                if (j == numVert - 1) {
                    mesh.getEdges().at(edgeIndex).setNextBElem(mesh.getCells().at(cellIndex).getBoundaryElementIndex(), cellIndex);
                }
                prevEdge = edgeIndex;
            }

        }
    }


    void loadCellTypes(std::istream& ist, MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;
        cellTypes.alocateData(mesh);
        for (IndexType i = 0; i < numCells; i++) {
            int vtkType = 0;
            ist >> vtkType;
            typename std::map<int, typename reader::type::ElementType>::iterator typeIt = TypeConversionTable.find(vtkType);
            if (typeIt != TypeConversionTable.end()){
                cellTypes.template getDataByPos<0>().at(i) = typeIt->second;
            } else {
                std::runtime_error("unsuported cell type");
            }
        }
    }


    void loadFromStream(std::istream& ist,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        ist.seekg(ist.beg);
        // Ignore first row "# vtk DataFile Version 2.0"
        ist.ignore(1024, '\n');
        // Ignore name of the data set
        ist.ignore(1024, '\n');
        // ASCII or BINARY
        std::string buf;
        std::getline(ist, buf);
        if (buf != "ASCII"){
            throw std::runtime_error("ASCII expected but got " + buf);
        }

        ist >> buf;
        if (buf != "DATASET"){
            throw std::runtime_error("the keyword DATASET expected");
        }

        ist >> buf;
        if (buf != "UNSTRUCTURED_GRID"){
            throw std::runtime_error("only unstructured grid is supported but got " + buf);
        }


        ist >> buf;
        if (buf == "POINTS") {
            loadPoints(ist, mesh);
        }

        ist >> buf;
        if (buf == "CELLS") {
            loadCells(ist, mesh);
        }

        ist >> buf;
        if (buf == "CELL_TYPES") {
            loadCellTypes(ist, mesh);
        }

    }

    MeshElements<3, IndexType, Real, Reserve...> loadFromStream(std::istream& ist){
        MeshElements<3, IndexType, Real, Reserve...> resultMesh;
        loadFromStream(ist, resultMesh);
        return resultMesh;
    }
};

#endif // VTKMESHREADER_H
