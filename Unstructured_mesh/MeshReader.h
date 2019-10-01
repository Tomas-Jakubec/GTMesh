#ifndef MESHREADER_H
#define MESHREADER_H


template<unsigned int MeshDimension, typename IndexType, typename Real>
class MeshReader{

};


template <typename IndexType, typename Real>
class MeshReader<2, IndexType, Real> {
public:
    enum ElementType{
        LINE = 1,
        TRIANGLE,
        QUAD,
        POLYGON
    };


};
#endif // MESHREADER_H
