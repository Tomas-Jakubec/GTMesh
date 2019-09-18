#ifndef COMPUTATIONALY_SIGNIFICANT_ELEMENT_H
#define COMPUTATIONALY_SIGNIFICANT_ELEMENT_H
#include "vertex.h"

template <unsigned int MeshDim, typename Real>
class ComputationallySignificantElement
{
protected:
    Vertex<MeshDim, Real> Center;
    int Flag;
public:
    ComputationallySignificantElement() {
        Center = {};
        Flag = int();
    }

    Vertex<MeshDim, Real>& GetCenter(){
        return Center;
    }

    const Vertex<MeshDim, Real>& GetCenter() const {
        return Center;
    }

    void SetCenter(const Vertex<MeshDim, Real>& v) {
        Center = v;
    }

    int& GetFlag() {
        return Flag;
    }

    const int& GetFlag() const {
        return Flag;
    }
};


#endif // COMPUTATIONALY_SIGNIFICANT_ELEMENT_H
