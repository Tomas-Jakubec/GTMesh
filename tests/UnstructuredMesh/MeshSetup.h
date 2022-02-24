// Test of several arithmetic operations
#include "GTMesh/UnstructuredMesh/UnstructuredMesh.h"


template <unsigned int ...Reserve>
void cube(UnstructuredMesh<3, size_t, double, Reserve...>& mesh3){

        mesh3.getVertices().push_back({0, {0,0,0}});
        mesh3.getVertices().push_back({1, {1,0,0}});
        mesh3.getVertices().push_back({2, {0,1,0}});
        mesh3.getVertices().push_back({3, {1,1,0}});
        mesh3.getVertices().push_back({4, {0,0,1}});
        mesh3.getVertices().push_back({5, {1,0,1}});
        mesh3.getVertices().push_back({6, {0,1,1}});
        mesh3.getVertices().push_back({7, {1,1,1}});

        mesh3.getEdges().push_back({0,0,1});
        mesh3.getEdges().push_back({1,0,2});
        mesh3.getEdges().push_back({2,1,3});
        mesh3.getEdges().push_back({3,2,3});
        mesh3.getEdges().push_back({4,0,4});
        mesh3.getEdges().push_back({5,1,5});
        mesh3.getEdges().push_back({6,3,7});
        mesh3.getEdges().push_back({7,2,6});
        mesh3.getEdges().push_back({8,4,5});
        mesh3.getEdges().push_back({9,5,7});
        mesh3.getEdges().push_back({10,6,7});
        mesh3.getEdges().push_back({11,4,6});

        mesh3.getFaces().push_back(0);
        mesh3.getFaces().at(0).getSubelements().addSubelement(0);
        mesh3.getFaces().at(0).getSubelements().addSubelement(1);
        mesh3.getFaces().at(0).getSubelements().addSubelement(2);
        mesh3.getFaces().at(0).getSubelements().addSubelement(3);
        mesh3.getFaces().push_back(1);
        mesh3.getFaces().at(1).getSubelements().addSubelement(0);
        mesh3.getFaces().at(1).getSubelements().addSubelement(5);
        mesh3.getFaces().at(1).getSubelements().addSubelement(8);
        mesh3.getFaces().at(1).getSubelements().addSubelement(4);
        mesh3.getFaces().push_back(2);
        mesh3.getFaces().at(2).getSubelements().addSubelement(1);
        mesh3.getFaces().at(2).getSubelements().addSubelement(4);
        mesh3.getFaces().at(2).getSubelements().addSubelement(11);
        mesh3.getFaces().at(2).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(3);
        mesh3.getFaces().at(3).getSubelements().addSubelement(3);
        mesh3.getFaces().at(3).getSubelements().addSubelement(6);
        mesh3.getFaces().at(3).getSubelements().addSubelement(10);
        mesh3.getFaces().at(3).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(4);
        mesh3.getFaces().at(4).getSubelements().addSubelement(2);
        mesh3.getFaces().at(4).getSubelements().addSubelement(6);
        mesh3.getFaces().at(4).getSubelements().addSubelement(9);
        mesh3.getFaces().at(4).getSubelements().addSubelement(5);
        mesh3.getFaces().push_back(5);
        mesh3.getFaces().at(5).getSubelements().addSubelement(8);
        mesh3.getFaces().at(5).getSubelements().addSubelement(9);
        mesh3.getFaces().at(5).getSubelements().addSubelement(10);
        mesh3.getFaces().at(5).getSubelements().addSubelement(11);


        mesh3.getFaces().at(0).setNextBElem(1,0);
        mesh3.getFaces().at(1).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(4,0);
        mesh3.getFaces().at(4).setNextBElem(5,0);
        mesh3.getFaces().at(5).setNextBElem(0,0);

        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(3);
        mesh3.updateSignature();

}
template <unsigned int ...Reserve>
void twoPrisms(UnstructuredMesh<3, size_t, double, Reserve...>& mesh3){
        using MeshType = UnstructuredMesh<3, size_t, double, Reserve...>;
        mesh3.getVertices().push_back({0, {0,0,0}});
        mesh3.getVertices().push_back({1, {1,0,0}});
        mesh3.getVertices().push_back({2, {0,1,0}});
        mesh3.getVertices().push_back({3, {1,1,0}});
        mesh3.getVertices().push_back({4, {0,0,1}});
        mesh3.getVertices().push_back({5, {1,0,1}});
        mesh3.getVertices().push_back({6, {0,1,1}});
        mesh3.getVertices().push_back({7, {1,1,1}});

        mesh3.getEdges().push_back({0,0,1});
        mesh3.getEdges().push_back({1,0,2});
        mesh3.getEdges().push_back({2,1,3});
        mesh3.getEdges().push_back({3,2,3});
        mesh3.getEdges().push_back({4,0,4});
        mesh3.getEdges().push_back({5,1,5});
        mesh3.getEdges().push_back({6,3,7});
        mesh3.getEdges().push_back({7,2,6});
        mesh3.getEdges().push_back({8,4,5});
        mesh3.getEdges().push_back({9,5,7});
        mesh3.getEdges().push_back({10,6,7});
        mesh3.getEdges().push_back({11,4,6});
        mesh3.getEdges().push_back({12,6,5});
        mesh3.getEdges().push_back({13,2,1});

    size_t index = 0;
        mesh3.getFaces().push_back(typename MeshType::Face(index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().push_back(typename MeshType::Face(++index));
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);


        mesh3.getFaces().at(0).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(8,0);
        mesh3.getFaces().at(8).setNextBElem(6,0);
        mesh3.getFaces().at(6).setNextBElem(0,0);
        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(0);

        mesh3.getFaces().at(1).setNextBElem(5,1);
        mesh3.getFaces().at(5).setNextBElem(4,1);
        mesh3.getFaces().at(4).setNextBElem(8,1);
        mesh3.getFaces().at(8).setNextBElem(7,1);
        mesh3.getFaces().at(7).setNextBElem(1,1);
        mesh3.getCells().push_back(1);
        mesh3.getCells().at(1).setBoundaryElementIndex(1);

        mesh3.updateSignature();

}
template <unsigned int ...Reserve>
void twoDeformedPrisms(UnstructuredMesh<3, size_t, double, Reserve...>& mesh3){

        mesh3.getVertices().push_back({0, {0,0,1}});
        mesh3.getVertices().push_back({1, {1,0,0}});
        mesh3.getVertices().push_back({2, {0,1,0}});
        mesh3.getVertices().push_back({3, {1,1,0}});
        mesh3.getVertices().push_back({4, {0,0,2}});
        mesh3.getVertices().push_back({5, {1,0,1}});
        mesh3.getVertices().push_back({6, {0,1,1}});
        mesh3.getVertices().push_back({7, {1,1,2}});

        mesh3.getEdges().push_back({0,0,1});
        mesh3.getEdges().push_back({1,0,2});
        mesh3.getEdges().push_back({2,1,3});
        mesh3.getEdges().push_back({3,2,3});
        mesh3.getEdges().push_back({4,0,4});
        mesh3.getEdges().push_back({5,1,5});
        mesh3.getEdges().push_back({6,3,7});
        mesh3.getEdges().push_back({7,2,6});
        mesh3.getEdges().push_back({8,4,5});
        mesh3.getEdges().push_back({9,5,7});
        mesh3.getEdges().push_back({10,6,7});
        mesh3.getEdges().push_back({11,4,6});
        mesh3.getEdges().push_back({12,6,5});
        mesh3.getEdges().push_back({13,2,1});

        size_t index = 0;
        mesh3.getFaces().push_back(index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(0);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(1);
        mesh3.getFaces().at(index).getSubelements().addSubelement(4);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(3);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(2);
        mesh3.getFaces().at(index).getSubelements().addSubelement(6);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(8);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(11);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(9);
        mesh3.getFaces().at(index).getSubelements().addSubelement(10);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().push_back(++index);
        mesh3.getFaces().at(index).getSubelements().addSubelement(12);
        mesh3.getFaces().at(index).getSubelements().addSubelement(7);
        mesh3.getFaces().at(index).getSubelements().addSubelement(13);
        mesh3.getFaces().at(index).getSubelements().addSubelement(5);


        mesh3.getFaces().at(0).setNextBElem(2,0);
        mesh3.getFaces().at(2).setNextBElem(3,0);
        mesh3.getFaces().at(3).setNextBElem(8,0);
        mesh3.getFaces().at(8).setNextBElem(6,0);
        mesh3.getFaces().at(6).setNextBElem(0,0);
        mesh3.getCells().push_back(0);
        mesh3.getCells().at(0).setBoundaryElementIndex(0);

        mesh3.getFaces().at(1).setNextBElem(5,1);
        mesh3.getFaces().at(5).setNextBElem(4,1);
        mesh3.getFaces().at(4).setNextBElem(8,1);
        mesh3.getFaces().at(8).setNextBElem(7,1);
        mesh3.getFaces().at(7).setNextBElem(1,1);
        mesh3.getCells().push_back(1);
        mesh3.getCells().at(1).setBoundaryElementIndex(1);

        mesh3.updateSignature();

}

template <unsigned int ...Reserve>
void square(UnstructuredMesh<2, size_t, double, Reserve...>& mesh){

    mesh.getVertices().resize(4);

    mesh.getVertices().at(0).setIndex(0);
    mesh.getVertices().at(0) = {0,0};
    mesh.getVertices().at(1).setIndex(1);
    mesh.getVertices().at(1) = {0,1};
    mesh.getVertices().at(2).setIndex(2);
    mesh.getVertices().at(2) = {1,0};
    mesh.getVertices().at(3).setIndex(3);
    mesh.getVertices().at(3) = {1,1};

    mesh.getFaces().resize(5);
    mesh.getFaces().at(0).vertexAIndex = 1;
    mesh.getFaces().at(0).vertexBIndex = 0;
    mesh.getFaces().at(1).vertexAIndex = 0;
    mesh.getFaces().at(1).vertexBIndex = 2;
    mesh.getFaces().at(2).vertexAIndex = 2;
    mesh.getFaces().at(2).vertexBIndex = 1;
    mesh.getFaces().at(3).vertexAIndex = 3;
    mesh.getFaces().at(3).vertexBIndex = 1;
    mesh.getFaces().at(4).vertexAIndex = 2;
    mesh.getFaces().at(4).vertexBIndex = 3;
    for(size_t i = 0; i < 5; i++)
        mesh.getFaces().at(i).setIndex(i);

    mesh.getFaces().at(0).setNextBElem(1,0);
    mesh.getFaces().at(1).setNextBElem(2,0);
    mesh.getFaces().at(2).setNextBElem(0,0);

    mesh.getFaces().at(2).setNextBElem(3,1);
    mesh.getFaces().at(3).setNextBElem(4,1);
    mesh.getFaces().at(4).setNextBElem(2,1);


    mesh.getCells().resize(2);
    mesh.getCells().at(0).setBoundaryElementIndex(0);
    mesh.getCells().at(0).setIndex(0);
    mesh.getCells().at(1).setBoundaryElementIndex(2);
    mesh.getCells().at(1).setIndex(1);
}


template <typename IndexType, typename Real, unsigned int ... Reserve>
void Pyramid4D(UnstructuredMesh<4, IndexType, Real, Reserve...>& mesh3){


    mesh3.getVertices().push_back({0, {0,0,0,0}});
    mesh3.getVertices().push_back({1, {1,0,0,0}});
    mesh3.getVertices().push_back({2, {0,1,0,0}});
    mesh3.getVertices().push_back({3, {0,0,1,0}});
    mesh3.getVertices().push_back({4, {0,0,0,1}});




    size_t index = 0;
    mesh3.getEdges().push_back({index++,0,1});
    mesh3.getEdges().push_back({index++,0,2});
    mesh3.getEdges().push_back({index++,1,2});
    mesh3.getEdges().push_back({index++,0,3});
    mesh3.getEdges().push_back({index++,1,3});
    mesh3.getEdges().push_back({index++,2,3});

    mesh3.getEdges().push_back({index++,0,4});//6
    mesh3.getEdges().push_back({index++,1,4});//7
    mesh3.getEdges().push_back({index++,2,4});//8
    mesh3.getEdges().push_back({index++,3,4});//9

    // Dim 2
    index = 0;
    mesh3.template getElements<2>().push_back(index++);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(0);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(1);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(2);

    mesh3.template getElements<2>().push_back(index++);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(0);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(3);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(4);

    mesh3.template getElements<2>().push_back(index++);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(2);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(4);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(5);

    mesh3.template getElements<2>().push_back(index++);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(1);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(3);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(5);
    // postup do 4D
    mesh3.template getElements<2>().push_back(index++);//4
    mesh3.template getElements<2>().back().getSubelements().addSubelement(0);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(6);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(7);

    mesh3.template getElements<2>().push_back(index++);//5
    mesh3.template getElements<2>().back().getSubelements().addSubelement(1);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(6);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(8);

    mesh3.template getElements<2>().push_back(index++);//6
    mesh3.template getElements<2>().back().getSubelements().addSubelement(2);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(7);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(8);

    mesh3.template getElements<2>().push_back(index++);//7
    mesh3.template getElements<2>().back().getSubelements().addSubelement(3);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(6);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(9);

    mesh3.template getElements<2>().push_back(index++);//8
    mesh3.template getElements<2>().back().getSubelements().addSubelement(4);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(7);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(9);

    mesh3.template getElements<2>().push_back(index++);//9
    mesh3.template getElements<2>().back().getSubelements().addSubelement(5);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(8);
    mesh3.template getElements<2>().back().getSubelements().addSubelement(9);



    // 3D
    index = 0;
    mesh3.template getElements<3>().push_back(index++);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(0);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(1);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(2);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(3);

    mesh3.template getElements<3>().push_back(index++);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(0);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(4);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(5);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(6);

    mesh3.template getElements<3>().push_back(index++);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(1);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(4);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(7);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(8);

    mesh3.template getElements<3>().push_back(index++);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(2);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(6);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(8);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(9);

    mesh3.template getElements<3>().push_back(index++);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(3);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(5);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(7);
    mesh3.template getElements<3>().back().getSubelements().addSubelement(9);
    // 4D
    mesh3.getCells().push_back(0);

    for (auto& face : mesh3.getFaces()) {
        face.setCellIndex(0);
        face.setNextBElem((face.getIndex() + 1)%5, 0);
    }
    mesh3.getCells()[0].setBoundaryElementIndex(0);
    mesh3.updateSignature();

}
