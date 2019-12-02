#ifndef GRAMMSCHMIDT_H
#define GRAMMSCHMIDT_H


#include "Vector.h"
#include <array>

/**
 * @brief GrammSchmidt
 * Gramm-Schmidt process without normalization. Can be used to
 * calculate volume in any dimension.
 * @param vectors [in / out] array containing the vectors to be processed.
 * After the process the input vectors are changed.
 */
template <unsigned int NumVec, unsigned int Dimension,typename IndexType, typename Real>
void grammSchmidt(std::array<Vector<Dimension, Real>, NumVec>& vectors){
    /*
     * Vector of inverse suquare of norm
     */
    std::array<Real, NumVec> invSumOfSquares;
    invSumOfSquares.at(0) = 1.0 / vectors.at(0).sumOfSquares();
    std::array<Real, NumVec> coef;

    for (IndexType i = 1; i < vectors.size(); i++) {
        /*
         * Coefitiens of scalar products.
         * The coefitients are computed in advance for better stability
         */
        for (IndexType j = 0; j < i; j++) {

            coef.at(j) = (vectors.at(i)*(vectors.at(j))) * invSumOfSquares.at(j);

        }
        for (IndexType j = 0; j < i; j++) {

            vectors.at(i) -= (vectors.at(j) * coef.at(j));
        }

        invSumOfSquares.at(i) = 1.0 / vectors.at(i).sumOfSquares();
    }
}


/**
 * @brief GrammSchmidt
 * Gramm-Schmidt process with normalization. Can be used to
 * calculate volume in any dimension or normal vector. The norms of orthogonalized
 * vectors are returned using norms array. This function is slightly more efficient
 * than the method without nomalization.
 * @param vectors [in / out] array containing the vectors to be processed.
 * After the process the vectors are changed.
 * @param norms [out] norms of the processed vectors.
 */
template <unsigned int NumVec, unsigned int Dimension,typename IndexType, typename Real>
void grammSchmidt(std::array<Vector<Dimension, Real>, NumVec>& vectors, std::array<Real, NumVec>& norms){
    /*
     * Vector of inverse suquare of norm
     */
    std::array<Real, NumVec> coef;

    norms.at(0) = vectors.at(0).normEukleid();
    vectors.at(0) /= norms.at(0);

    for (IndexType i = 1; i < vectors.size(); i++) {
        /*
         * Coefitiens of scalar products.
         * The coefitients are computed in advance for better stability
         */
        for (IndexType j = 0; j < i; j++) {

            coef.at(j) = (vectors.at(i)*(vectors.at(j)));

        }
        for (IndexType j = 0; j < i; j++) {

            vectors.at(i) -= (vectors.at(j) * coef.at(j));

        }

        norms.at(i) = vectors.at(i).normEukleid();
        vectors.at(i) /= norms.at(i);
    }
}

#endif // GRAMMSCHMIDT_H
