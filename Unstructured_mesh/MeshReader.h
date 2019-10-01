#ifndef MESHREADER_H
#define MESHREADER_H


template<unsigned int MeshDimension>
class MeshReader{

};


template <>
class MeshReader<2> {
public:
    enum ElementType{
        LINE = 1,
        TRIANGLE,
        QUAD,
        POLYGON
    };


};
#endif // MESHREADER_H
