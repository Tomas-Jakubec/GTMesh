#ifndef GRAMMSCHMIDT_H
#define GRAMMSCHMIDT_H


#include "../../NumericStaticArray/Vertex.h"
#include <vector>

/**
 * @brief GrammSchmidt
 * Gramm-Schmidt process without normalization. Can be used to
 * calculate volume in any dimension.
 * @param vectors [in / out] vector containing the vectors and
 * after the process the vectors are changed.
 */
template <unsigned int Dimension,typename IndexType, typename Real>
void GrammSchmidt(std::vector<Vertex<Dimension, Real>>& vectors){
    /*
     * Vector of inverse suquare of norm
     */
    std::vector<Real> invSumOfSquares(vectors.size());
    invSumOfSquares.at(0) = 1.0 / vectors.at(0).sumOfSquares();

    for (IndexType i = 1; i < vectors.size(); i++) {
        /*
         * Coefitiens of scalar products.
         * The coefitients are computed in advance for better stability
         */
        std::vector<Real> coef(i);
        for (IndexType j = 0; j < i; j++) {

            coef.at(j) = (vectors.at(i)*(vectors.at(j))) * invSumOfSquares.at(j);

        }
        for (IndexType j = 0; j < i; j++) {

            vectors.at(i) -= (vectors.at(j) * coef.at(j));
        }

        invSumOfSquares.at(i) = 1.0 / vectors.at(i).sumOfSquares();
    }
}

#endif // GRAMMSCHMIDT_H
