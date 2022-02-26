#ifndef MESHELEMENTS_H
#define MESHELEMENTS_H

#include "MeshElement.h"
#include "ElementIndex.hpp"

template <unsigned int Dimension, typename IndexType, typename Real, unsigned int ...Reserve>
struct MeshElements{
private:

    template<unsigned int dim, typename Void = void>
    struct _Reserve{

        static unsigned int constexpr value = std::get<Dimension - dim - 1>(std::array<unsigned int, sizeof... (Reserve)>{Reserve...});

    };

    template<unsigned int dim>
    struct _Reserve<dim, typename std::enable_if<dim == Dimension || dim == 1 || dim == 0 || (Dimension - dim > sizeof...(Reserve))>::type>{

        static unsigned int constexpr value = 0;
    };

public:

    template<unsigned int dim>
    static unsigned int constexpr reserve() {
        return _Reserve<dim>::value;
    }

    using Vertex = MeshElement<Dimension, 0, IndexType, Real, 0>;
    using Edge = MeshElement<Dimension, 1, IndexType, Real, 0>;
    using Face = MeshElement<Dimension, Dimension - 1, IndexType, Real, _Reserve<Dimension - 1>::value>;
    using Cell = MeshElement<Dimension, Dimension, IndexType, Real, 0>;

    template<unsigned int ElementDimension>
    using ElementType = MeshElement<Dimension, ElementDimension, IndexType, Real, _Reserve<ElementDimension>::value>;


    static unsigned int constexpr meshDimension() {
        return Dimension;
    }

private:
    template <unsigned int ElemDim = Dimension, typename Dummy = void>
    struct _MeshElements : public _MeshElements<ElemDim - 1, Dummy>{
        std::vector<ElementType<ElemDim>> elements;
    };

    template <typename Dummy>
    struct _MeshElements<0, Dummy>{
        std::vector<ElementType<0>> elements;
    };


private:
    _MeshElements<Dimension> innerElements;
    std::vector<Cell> BoundaryCells;
    //_MeshElements<Dimension> boundaryElements;

    /**
     * @brief Hash signature of the mash elements.
     * Use to detect changes in mesh.
     */
    size_t meshSignature;


public:
    template<unsigned int dim>
    std::vector<ElementType<dim>>& getElements(){
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return innerElements._MeshElements<dim>::elements;
    }


    template<unsigned int dim>
    const std::vector<ElementType<dim>>&  getElements() const {
        static_assert (Dimension >= dim, "In GetElements template parameter dim must be less or equal to Dimension.");
        return innerElements._MeshElements<dim>::elements;
    }

    template<unsigned int dim>
    ElementType<dim>&  getElement(const ElementIndex<dim, IndexType>& elementIndex) {
        return getElements<dim>()[elementIndex.index];
    }

    template<unsigned int dim>
    const ElementType<dim>&  getElement(const ElementIndex<dim, IndexType>& elementIndex) const {
        return getElements<dim>()[elementIndex];
    }

    std::vector<Vertex>& getVertices(){
        return getElements<0>();
    }


    std::vector<Edge>& getEdges(){
        return getElements<1>();
    }

    std::vector<Face>& getFaces(){
        return getElements<Dimension - 1>();
    }

    std::vector<Cell>& getCells(){
        return getElements<Dimension>();
    }

    std::vector<Cell>& getBoundaryCells() {
        return BoundaryCells;
    }
/*
 * Constant version of getters
 */
    const std::vector<Vertex>& getVertices() const {
        return getElements<0>();
    }


    const std::vector<Edge>& getEdges() const {
        return getElements<1>();
    }

    const std::vector<Face>& getFaces() const {
        return getElements<Dimension - 1>();
    }

    const std::vector<Cell>& getCells() const {
        return getElements<Dimension>();
    }

    const std::vector<Cell>& getBoundaryCells() const {
        return BoundaryCells;
    }

private:
    template<unsigned int Dim, typename Dummy = void>
    struct _ClearMesh {
        static void clear(MeshElements<Dimension, IndexType, Real, Reserve...> &mesh){
            mesh.template getElements<Dim>().clear();
            _ClearMesh<Dim - 1>::clear(mesh);
        }
    };


    template<typename Dummy>
    struct _ClearMesh<0, Dummy> {
        static void clear(MeshElements<Dimension, IndexType, Real, Reserve...> &mesh){
            mesh.template getElements<0>().clear();
        }
    };

public:
    /**
     * @brief clear<HR>
     * Sets the size of all vectors in the mesh to 0 by calling clear on them.
     */
    void clear() {
        _ClearMesh<Dimension>::clear(*this);
    }

    void appendBoundaryCell(IndexType cellIndex, IndexType faceIndex){
        Cell c;
        c.setIndex(cellIndex);
        c.setBoundaryElementIndex(faceIndex);
        BoundaryCells.push_back(c);
    }

    void setupBoundaryCells(){
        for (Face& face : getFaces()){
            if (isInvalidIndex(face.getCellLeftIndex())){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.setCellLeftIndex(cellIndex);
                appendBoundaryCell(cellIndex, face.getIndex());
            }
            if (isInvalidIndex(face.getCellRightIndex())){
                IndexType cellIndex = BoundaryCells.size() | BOUNDARY_INDEX(IndexType);
                face.setCellRightIndex(cellIndex);
                appendBoundaryCell(cellIndex, face.getIndex());
            }
        }
        BoundaryCells.shrink_to_fit();
    }


    void setupBoundaryCellsCenters() {
        for(Cell& cell : BoundaryCells){
            cell.setCenter(getFaces().at(cell.getBoundaryElementIndex()).getCenter());
        }
    }




    /*
     * Signature computation
     */
private:

    template<unsigned int Dim = Dimension, typename Void = void>
    struct HashOfMeshElements{


        static size_t hash(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh)
        {
            std::hash<IndexType> indexHasher;

            // Hash of generic element eg. face
            size_t elemHash = mesh.getElements<Dim>().size();
            for(auto& element : mesh.getElements<Dim>()) {
                for (auto& subElement : element.getSubelements()) {
                    elemHash ^= indexHasher(subElement);
                }
            }
            return elemHash ^ HashOfMeshElements<Dim - 1>::hash(mesh);
        }


    };

    template<unsigned int Dim>
    struct HashOfMeshElements <Dim, typename std::enable_if<(MeshElements<Dimension, IndexType, Real, Reserve...>::template reserve<Dim>() > 0), size_t>::type>{
        static size_t hash(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh)
        {
#if __cplusplus <= 201702L // standard c++14 and older
            std::hash<std::string> hasher;
            // Use string as a byte container representing the array
            std::string tmpString(reinterpret_cast<char*>(mesh.template getElements<Dim>().data()), mesh.template getElements<Dim>().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::ElementType<Dim>));
            size_t elemHash = hasher(tmpString);
#else // standard c++17 and later
            std::hash<std::string_view> hasher;
            // Use string as a byte container representing the array
            std::string_view vectorView(reinterpret_cast<char*>(mesh.template getElements<Dim>().data()), mesh.template getElements<Dim>().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::ElementType<Dim>));
            size_t elemHash = hasher(vectorView);
#endif
            return elemHash ^ HashOfMeshElements<Dim - 1>::hash(mesh);
        }
    };


    /**
     * @brief The hashOfMeshElements<Dimension, Dummy> struct hashes the cell information.
     */
    template<typename Dummy>
    struct HashOfMeshElements<Dimension, Dummy>{
        static size_t hash(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

            // Hash of cells
#if __cplusplus <= 201702L // standard c++14 and older
            std::hash<std::string> hasher;
            // Use string as a byte container representing the array
            std::string tmpString(reinterpret_cast<char*>(mesh.getCells().data()), mesh.getCells().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::Cell));

            size_t cHash = hasher(tmpString);
#else // standard c++17 and later
            std::hash<std::string_view> hasher;
            // Use string as a byte container representing the array
            std::string_view vectorView(reinterpret_cast<char*>(mesh.getCells().data()), mesh.getCells().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::Cell));

            size_t cHash = hasher(vectorView);
#endif

            return cHash ^ HashOfMeshElements<Dimension -1>::hash(mesh);
        }
    };

    template<typename Dummy>
    struct HashOfMeshElements<1, Dummy>{
        static size_t hash(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

            // Hash of edges
#if __cplusplus <= 201702L // standard c++14 and older
            std::hash<std::string> hasher;
            // Use string as a byte container representing the array
            std::string tmpString(reinterpret_cast<char*>(mesh.getEdges().data()), mesh.getEdges().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::Edge));

            size_t eHash = hasher(tmpString);
#else // standard c++17 and later
            std::hash<std::string_view> hasher;
            // Use string as a byte container representing the array
            std::string_view vectorView(reinterpret_cast<char*>(mesh.getEdges().data()), mesh.getEdges().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::Edge));

            size_t eHash = hasher(vectorView);
#endif
            return eHash ^ HashOfMeshElements<0>::hash(mesh);
        }
    };


    template<typename Dummy>
    struct HashOfMeshElements<0, Dummy>{
        static size_t hash(MeshElements<Dimension, IndexType, Real, Reserve...>& mesh){

            // Hash of vertices
#if __cplusplus <= 201702L // standard c++14 and older
            std::hash<std::string> hasher;

            // Use string as a byte container representing the array
            std::string tmpString(reinterpret_cast<char*>(mesh.getVertices().data()), mesh.getVertices().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::Vertex));

            size_t vHash = hasher(tmpString);
#else // standard c++17 and later
            std::hash<std::string_view> hasher;

            // Use string as a byte container representing the array
            std::string_view tmpString(reinterpret_cast<char*>(mesh.getVertices().data()), mesh.getVertices().size() * sizeof (MeshElements<Dimension, IndexType, Real, Reserve...>::Vertex));

            size_t vHash = hasher(tmpString);
#endif
            return vHash;
        }
    };

public:

    size_t updateSignature() {

        meshSignature = HashOfMeshElements<Dimension>::hash(*this);

        return meshSignature;

    }

    size_t getSignature() const {
        return meshSignature;
    }
};


#endif // MESHELEMENTS_H
