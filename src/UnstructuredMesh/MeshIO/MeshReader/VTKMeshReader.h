#ifndef VTKMESHREADER_H
#define VTKMESHREADER_H

#include "MeshReader.h"
#include "../../MeshDataContainer/MeshDataContainer.h"
#include "../../MeshElements/MeshElements.h"
#include <istream>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>

/**
 *
 */
template<unsigned int MeshDimension>
class VTKMeshReader : public MeshReader<MeshDimension>{
public:
    VTKMeshReader() = default;

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    VTKMeshReader(const MeshElements<MeshDimension, IndexType, Real, Reserve...>&){}
};


/**
 *
 */
template<>
class VTKMeshReader<2> : public MeshReader<2>{
    using reader = MeshReader<2>;
    std::map<int, typename reader::type::ElementType> TypeConversionTable{
        {3, reader::type::ElementType::LINE},
        {5, reader::type::ElementType::TRIANGLE},
        {8, reader::type::ElementType::QUAD},
        {9, reader::type::ElementType::QUAD},
        {7, reader::type::ElementType::POLYGON},
    };


    MeshDataContainer<typename reader::type::ElementType, 2> cellTypes;
    // file indexing
    //

public:
    VTKMeshReader() = default;

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    VTKMeshReader(const MeshElements<2, IndexType, Real, Reserve...>&){}

    virtual ~VTKMeshReader() = default;

    virtual MeshDataContainer<typename reader::type::ElementType, 2> getCellTypes() const {
        return cellTypes;
    }

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadVertices(std::istream& ist, MeshElements<2, IndexType, Real, Reserve...>& mesh){
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

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadCells(std::istream& ist, MeshElements<2, IndexType, Real, Reserve...>& mesh){
        std::unordered_map<std::string, IndexType> edges;

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
                    edges[edgeKey] = edgeIndex;
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

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadCellTypes(std::istream& ist, MeshElements<2, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;
        cellTypes.allocateData(mesh);
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

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadFromStream(std::istream& ist,MeshElements<2, IndexType, Real, Reserve...>& mesh){
        ist.seekg(ist.beg);
        // Ignore first row "# vtk DataFile Version 2.0"
        ist.ignore(1024, '\n');
        // Ignore name of the data set
        ist.ignore(1024, '\n');
        // ASCII or BINARY
        std::string buf;

        ist >> buf;
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

        mesh.clear();

        ist >> buf;
        if (buf == "POINTS") {
            loadVertices(ist, mesh);
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

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    MeshElements<2, IndexType, Real, Reserve...> loadFromStream(std::istream& ist){
        MeshElements<2, IndexType, Real, Reserve...> resultMesh;
        loadFromStream(ist, resultMesh);
        return resultMesh;
    }
};

/**
 *
 */
template<>
class VTKMeshReader<3> : public MeshReader<3>{
    using reader = MeshReader<3>;
    std::map<int, typename reader::type::ElementType> TypeConversionTable{
        {10, reader::type::ElementType::TETRA},
        {11, reader::type::ElementType::HEXAHEDRON},
        {12, reader::type::ElementType::HEXAHEDRON},
        {13, reader::type::ElementType::WEDGE},
        {14, reader::type::ElementType::PYRAMID},
    };
    const std::map<
        const int,
        const std::pair<
            const std::vector<std::array<int,2>>,
            const std::vector<std::vector<int>
            >
        >
    > TypeEdgesFaces{
        {4, {// tetrahedron
                {// edges (first)
                    {0,1},{1,2},{2,0},{0,3},{1,3},{2,3}
                },{//faces (second)
                    {0,1,2}, {0,4,3}, {1,5,4}, {3,5,2}
                }
            }
        },
        {8, {// hexahedron
                {// edges (first)
                    {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},{2,6},{3,7},{4,5},{5,6},{6,7},{7,4}
                },{//faces (second)
                    {0,1,2,3},
                    {4,0,5,8},
                    {5,1,6,9},
                    {6,2,7,10},
                    {7,3,5,11},
                    {8,9,10,11}
                }
            }
        },
        {6, {// wedge
                {// edges (first)
                    {0,1},
                    {1,2},
                    {2,0},
                    {0,3},
                    {1,4},
                    {2,5},
                    {3,5},
                    {5,4},
                    {4,3}
                },{//faces (second)
                    {0,2,1},
                    {0,4,8,3},
                    {1,4,7,5},
                    {3,5,6,2},
                    {6,7,8}
                }
            }
        },
        {5, {// pyramid
                {// edges (first)
                    {0,1},//0
                    {1,2},//1
                    {2,3},//2
                    {3,0},//3
                    {0,4},//4
                    {1,4},//5
                    {2,4},//6
                    {3,4} //7
                },{//faces (second)
                    {0,1,2,3},
                    {0,5,4},
                    {1,6,5},
                    {2,7,6},
                    {3,4,7}
                }
            }
        },
    };

    MeshDataContainer<typename reader::type::ElementType, 3> cellTypes;
    // file indexing
    //

    //
    //MeshDataContainer<IndexType>
public:
    VTKMeshReader() = default;

    virtual ~VTKMeshReader() = default;

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    VTKMeshReader(const MeshElements<3, IndexType, Real, Reserve...>&){}

    virtual MeshDataContainer<typename reader::type::ElementType, 3> getCellTypes() const {
        return cellTypes;
    }

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadVertices(std::istream& ist, MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numPoints;
        ist >> numPoints;
        mesh.getVertices().resize(numPoints);
        ist.ignore(20, '\n');
        for (IndexType vertIndex = 0; vertIndex < numPoints; vertIndex++) {
            mesh.getVertices().at(vertIndex).setIndex(vertIndex);
            ist >> mesh.getVertices().at(vertIndex)[0];
            ist >> mesh.getVertices().at(vertIndex)[1];
            ist >> mesh.getVertices().at(vertIndex)[2];
        }
    }

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadCells(std::istream& ist, MeshElements<3, IndexType, Real, Reserve...>& mesh){

        std::unordered_map<std::string, IndexType> edges;
        std::unordered_map<std::string, IndexType> faces;

        IndexType numCells;
        ist >> numCells;
        mesh.getCells().resize(numCells);
        // Skip the total number of numbers
        ist.ignore(50, '\n');
        for (IndexType cellIndex = 0; cellIndex < numCells; cellIndex++) {
            mesh.getCells().at(cellIndex).setIndex(cellIndex);
            IndexType numVert;
            ist >> numVert;
            // load vertices of the cell
            std::vector<IndexType> vertices(numVert);
            for(IndexType j = 0; j < numVert; j++){
                ist >> vertices.at(j);
            }

            // construct an element
            // obtain constructing order of edges and faces
            const std::vector<std::array<int,2>>& edgeOrder = TypeEdgesFaces.at(numVert).first;
            const std::vector<std::vector<int>>& faceOrder = TypeEdgesFaces.at(numVert).second;

            std::vector<IndexType> edgeIndexes;
            // construct edges first
            for (const std::array<int, 2>& e : edgeOrder) {

                IndexType iA = vertices.at(e[0]), iB = vertices.at(e[1]);
                std::string edgeKey = iA < iB ? std::to_string(iA) +";"+ std::to_string(iB) : std::to_string(iB) +";"+ std::to_string(iA);
                typename std::unordered_map<std::string, IndexType>::iterator edgeIt = edges.find(edgeKey);

                if (edgeIt == edges.end()) {
                    IndexType edgeIndex = IndexType();
                    edgeIndex = mesh.getEdges().size();
                    mesh.getEdges().push_back({});
                    mesh.getEdges().at(edgeIndex).setVertexAIndex(iA);
                    mesh.getEdges().at(edgeIndex).setVertexBIndex(iB);
                    mesh.getEdges().at(edgeIndex).setIndex(edgeIndex);
                    edgeIndexes.push_back(edgeIndex);
                    edges[edgeKey] = edgeIndex;
                } else {
                    edgeIndexes.push_back(edgeIt->second);
                }
            }

            IndexType prevFaceIndex = INVALID_INDEX(IndexType);
            for (IndexType fi = 0; fi < faceOrder.size(); fi++) {

                const std::vector<int>& f = faceOrder.at(fi);

                std::vector<IndexType> faceEdges;
                for (const int& index : f) {
                    faceEdges.push_back(edgeIndexes.at(index));
                }
                std::sort(faceEdges.begin(), faceEdges.end());

                std::string faceKey = "";
                for (IndexType& eI : faceEdges) {
                    faceKey += std::to_string(eI) + ";";
                }

                typename std::unordered_map<std::string, IndexType>::iterator faceIt = faces.find(faceKey);

                IndexType faceIndex;
                if (faceIt == faces.end()) {
                    faceIndex = mesh.getFaces().size();
                    mesh.getFaces().push_back({});
                    for (const int& index : f) {
                        mesh.getFaces().at(faceIndex).getSubelements().addSubelement(edgeIndexes.at(index));
                    }
                    mesh.getFaces().at(faceIndex).setCellLeftIndex(cellIndex);
                    mesh.getFaces().at(faceIndex).setIndex(faceIndex);
                } else {
                    faceIndex = faceIt->second;
                    mesh.getFaces().at(faceIndex).setCellRightIndex(cellIndex);
                }

                if (prevFaceIndex != INVALID_INDEX(IndexType)) {
                    mesh.getFaces().at(prevFaceIndex).setNextBElem(faceIndex, cellIndex);
                }
                if (fi == 0) {
                    mesh.getCells().at(cellIndex).setBoundaryElementIndex(faceIndex);
                }
                if (fi == faceOrder.size() - 1) {
                    mesh.getFaces().at(faceIndex).setNextBElem(mesh.getCells().at(cellIndex).getBoundaryElementIndex(), cellIndex);
                }
                faces[faceKey] = faceIndex;
                prevFaceIndex = faceIndex;
            }
        }
    }

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadCellTypes(std::istream& ist, MeshElements<3, IndexType, Real, Reserve...>& mesh){
        IndexType numCells;
        ist >> numCells;
        cellTypes.allocateData(mesh);
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

    template<typename IndexType, typename Real, unsigned int ...Reserve>
    void loadFromStream(std::istream& ist,MeshElements<3, IndexType, Real, Reserve...>& mesh){
        ist.seekg(ist.beg);
        // Ignore first row "# vtk DataFile Version 2.0"
        ist.ignore(1024, '\n');
        // Ignore name of the data set
        ist.ignore(1024, '\n');
        // ASCII or BINARY
        std::string buf;
        ist >> buf;

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

        mesh.clear();

        ist >> buf;
        if (buf == "POINTS") {
            loadVertices(ist, mesh);
        }

        ist >> buf;
        if (buf == "CELLS") {
            loadCells(ist, mesh);
        }

        ist >> buf;
        if (buf == "CELL_TYPES") {
            loadCellTypes(ist, mesh);
        }

        mesh.updateSignature();
    }


    template<typename IndexType, typename Real, unsigned int ...Reserve>
    MeshElements<3, IndexType, Real, Reserve...> loadFromStream(std::istream& ist){
        MeshElements<3, IndexType, Real, Reserve...> resultMesh;
        loadFromStream(ist, resultMesh);
        return resultMesh;
    }
};

#endif // VTKMESHREADER_H
