/**
\file AffineState.h
\brief All functions needed for affine state-based simulation described in arXiv:2109.08629
*/

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>

#include <qpp/qpp.h>

using mat_u_t = Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic>;
using vec_u_t = Eigen::Vector<unsigned, Eigen::Dynamic>;
using block_t = Eigen::Block<mat_u_t>;

namespace stab {
    class AffineState {
    public:

        explicit AffineState(int n);

        AffineState(const AffineState &other);

        void CZ(int a, int b); /*! \brief Controlled-\f$ Z \f$ gate*/
        void CX(int a, int b); /*! \brief Controlled-\f$ X \f$ gate
                               * \param a Control qubit
                               * \param b Target qubit
                               */
        void SWAP(int a, int b); /*! \brief Swap gate*/
        void S(int j); /*! \brief The gate \f$S = |0\rangle\langle0| + i|1\rangle\langle1|\f$ */
        void SDG(int j); /*! \brief The gate \f$S^\dagger = |0\rangle\langle0| - i|1\rangle\langle1|\f$*/
        void H(int a);   /*! \brief Hadamard gate */
        void X(int j); /*! \brief Pauli \f$ X \f$ gate*/
        void Y(int j); /*! \brief Pauli \f$ Y \f$ gate*/
        void Z(int j); /*! \brief Pauli \f$ Z \f$ gate*/
        int MeasureZ(int j, bool postselect = false,
                     int postselected_outcome = 0);  /*! \brief Measure a single qubit
                                                     * \param j Index of qubit to be measured
                                                     * \param postselect Optional flag indicating whether postselection should be used
                                                     * \param postselected_outcome Desired measurement outcome if postselection is used. Throws error if postselection is chosen but is impossible. 
                                                     * \return Measurement outcome (0 or 1)
                                                     */
        void Reset(int j); /*! \brief Resets qubit to \f$|0\rangle\f$*/
        std::vector<int> MeasureAll() const; /*! \brief Measure all qubits
                                             * \return Vector of measurement outcomes
                                             */
        std::map<std::vector<int>, int> Sample(int nreps) const; /*! \brief Sample the state repeatedly
                                                                 * \param nreps Number of samples
                                                                 * \return Map whose keys are vectors of measurement outcomes and whose values are number of observations with each result
                                                                 */
        Eigen::VectorXcd to_vec() const; /*!\brief Construct full statevector*/

        // Get member variables:
        int n() const; /*! \brief Get number of qubits*/

        int phase() const;  /*! \brief Get value \f$ p \f$ such that global phase is \f$ e^{ip\pi/4}\f$ */

        mat_u_t Q() const;  /*! \brief Get \f$ Q \f$ matrix containing phases */

        mat_u_t A() const;  /*! \brief Get \f$ A \f$ generating the affine space*/

        vec_u_t b() const; /*! \brief Get offset vector \f$ b \f$*/

        std::vector<int> pivots() const; /*! \brief Get the principal index map. Entry \f$ i \f$ corresponds to the value \f$ p(i) \f$.*/

        int r() const; /*! \brief Get rank \f$ r \f$ of generating matrix*/

        friend std::ostream &operator<<(std::ostream &out, AffineState const &psi);

    private:
        int n_, phase_, r_;     // number of qubits
        mat_u_t Q_, A_;
        vec_u_t b_;
        std::vector<int> pivots_;

        // Subroutines
        std::vector<int> A_col_nonzeros(int row);

        std::vector<int> A_row_nonzeros(int col);

        std::vector<int> Q_nonzeros();

        void FixFinalBit(int z);

        void ReduceGramRowCol(int c);

        void ReindexSubtColumn(int k, int c, std::vector<int> col_c_nonzeros);

        void ReindexSwapColumns(int k);

        void MakePrincipal(int c, int j);

        bool ReselectPrincipalRow(int j, int c);

        void ZeroColumnElim(int c);

        void S_or_SDG(int j, bool dg);

        int piv_col(int row_number);

    }; // class AffineState

    template<typename Derived>
    using expr_t =
            typename std::decay<decltype(std::declval<Derived>().eval())>::type;

    // Small helper functions
    auto mod_2 = []() {
        // returns a unary lambda
        return [](unsigned x) { return x & 1; };
    };

    auto inline m2 = mod_2();

    template<typename Derived>
    [[nodiscard]] expr_t<Derived> ReduceMod2(const Eigen::MatrixBase<Derived> &A) {
        return A.unaryExpr(m2);
    }

} // namespace stab
