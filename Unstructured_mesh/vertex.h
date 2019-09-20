#ifndef VERTEX_H
#define VERTEX_H
#include <math.h>
#include <iostream>
#include <initializer_list>
#include "inline_array_operations.h"


template <unsigned int Dim, typename Real = double>
class Vertex {
    /**
     * @brief Coordinates
     */
    Real Coordinates[Dim] = {};

public:
    Vertex(){}
    Vertex(std::initializer_list<Real> l){
        *this = l;
    }

    Vertex<Dim, Real>& operator =(std::initializer_list<Real> l);

    void SetCoordinate(Real coord, unsigned int pos){
        Coordinates[pos] = coord;
    }

    Real& operator[](unsigned int pos){
        return Coordinates[pos];
    }

    const Real& operator[](unsigned int pos) const {
        return Coordinates[pos];
    }

    Real NormEukleid();

    inline Real SumOfSquares() {
        return  inlineScalarProduct<Dim, Real>::computation(Coordinates, Coordinates);
    }

    Vertex<Dim, Real> operator-(const Vertex<Dim, Real>&) const;

    Vertex<Dim, Real> operator+(const Vertex<Dim, Real>&) const;

    Vertex<Dim, Real> operator*(const Real&) const;

    Vertex<Dim, Real> operator/(const Real&) const;

    /**
     * @brief operator *
     *
     * Scalar product of two vectors
     * @param v
     * @return
     */
    Real operator*(const Vertex<Dim, Real>& v);


    Vertex<Dim, Real>& operator+=(const Vertex<Dim, Real>&);
    Vertex<Dim, Real>& operator-=(const Vertex<Dim, Real>&);
    Vertex<Dim, Real>& operator*=(const Real&);
    Vertex<Dim, Real>& operator/=(const Real&);

    bool operator==(const Vertex<Dim, Real>&) const;
    bool operator!=(const Vertex<Dim, Real>&) const;
};


template<unsigned int Dim, typename Real>
Vertex<Dim, Real>& Vertex<Dim, Real>::operator =(std::initializer_list<Real> l){
    unsigned int i = 0;

    for(Real x : l){
        if (i < Dim){
            Coordinates[i] = x;
        }else{
            break;
        }
        i++;
    }
    if (i < Dim){
        for (; i < Dim; i++) {
            Coordinates[i] = Real();
        }
    }
    return *this;
}
/*
** Calculates the Eucleid norm of the point
*/
template <unsigned int Dim, typename Real>
Real Vertex<Dim, Real>::NormEukleid(){
    return sqrt(SumOfSquares());
}




/*
** Overloaded operators
*/
//subtracting two points
template <unsigned int Dim, typename Real>
Vertex<Dim, Real> Vertex<Dim, Real>::operator -(const Vertex<Dim, Real>& v) const {
    Vertex<Dim, Real> res;
    inlineSubtraction<Dim, Real>::computation(res.Coordinates, this->Coordinates, v.Coordinates);
    return res;
}

//addition of two points
template <unsigned int Dim, typename Real>
Vertex<Dim, Real> Vertex<Dim, Real>::operator +(const Vertex<Dim, Real>& v) const {
    Vertex<Dim, Real> res;
    inlineAddition<Dim, Real>::computation(res.Coordinates, this->Coordinates, v.Coordinates);
    return res;
}

//multiplying coordinates with real number
template <unsigned int Dim, typename Real>
Vertex<Dim, Real> Vertex<Dim, Real>::operator *(const Real& x) const {
    Vertex<Dim, Real> res;
    inlineMultiplication<Dim, Real>::computation(res.Coordinates, this->Coordinates, x);
    return res;
}


//division
template <unsigned int Dim, typename Real>
Vertex<Dim, Real> Vertex<Dim, Real>::operator /(const Real& x) const {
    return this->operator*(Real(1.0)/x);
}


template<unsigned int Dim, typename Real>
Real Vertex<Dim, Real>::operator*(const Vertex<Dim, Real> &v)
{
    return inlineScalarProduct<Dim, Real>::computation(Coordinates, v.Coordinates);
}



// Adds value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vertex<Dim, Real>& Vertex<Dim, Real>::operator +=(const Vertex<Dim, Real>& v){
    inlineAddition<Dim, Real>::computation(Coordinates, v.Coordinates);
    return *this;
}

// Subtracts value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vertex<Dim, Real>& Vertex<Dim, Real>::operator -=(const Vertex<Dim, Real>& v){
    inlineSubtraction<Dim, Real>::computation(Coordinates, v.Coordinates);
    return *this;
}


// Adds value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vertex<Dim, Real>& Vertex<Dim, Real>::operator *=(const Real& x){
    inlineMultiplication<Dim, Real>::computation(Coordinates, x);
    return *this;
}

// Subtracts value of coordinates of another Point
template <unsigned int Dim, typename Real>
Vertex<Dim, Real>& Vertex<Dim, Real>::operator /=(const Real& x){
    this->operator*=(Real(1.0)/x);
    return *this;
}

// Compares two points wether they are the same
template <unsigned int Dim, typename Real>
bool Vertex<Dim, Real>::operator ==(const Vertex<Dim, Real>& v) const {
    for(unsigned int i = 0; i < Dim; i++) {
        if(this->operator[](i) != v[i]){
            return false;
        }
    }
    return true;
}

// Compares two points wether they are not the same
template <unsigned int Dim, typename Real>
bool Vertex<Dim, Real>::operator !=(const Vertex<Dim, Real>& v) const {
    return !(*this == v);
}

template <unsigned int Dim, typename Real = double>
std::ostream& operator <<(std::ostream& ost, const Vertex<Dim,Real>& v) {
    for (unsigned int i = 0; i < Dim; i++) {
        ost << v[i] << ' ';
    }

    return ost;
}


template <unsigned int Dim, typename Real = double>
std::istream& operator >>(std::istream& ist, Vertex<Dim,Real>& v) {
    for (unsigned int i = 0; i < Dim; i++) {
        ist >> v[i];
    }
    return ist;
}


#endif // VERTEX_H
