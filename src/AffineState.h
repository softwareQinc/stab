/**
\file AffineState.h
\brief All functions needed for the Clifford simulation method described in Niel
de Beaudrap and Steven Herbert's paper "Fast Stabiliser Simulation with
Quadratic Form Expansions"
(https://quantum-journal.org/papers/q-2022-09-15-803/). Most functions are taken
directly from the pseudocode described in that paper. We remark that another
implementation of this algorithm is given in https://github.com/CQCL/simplex,
and that this was the source of the idea to, instead of iterating over an entire
row/column, first search for the nonzero entries and then iterate only over
them. This saves some time when doing nested for loops, and might give a small
speedup for single for loops when A and Q are sparse. \author softwareQ
*/

#ifndef STAB_AFFINE_STATE_H_
#define STAB_AFFINE_STATE_H_

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>

#ifdef USE_QPP
#include <qpp/qpp.h>
#endif // USE_QPP

using mat_u_t = Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic>;
using vec_u_t = Eigen::Vector<unsigned, Eigen::Dynamic>;
using block_t = Eigen::Block<mat_u_t>;

namespace stab {
class AffineState {
  public:
    /*! \brief Initialize a quantum state given by \f$|0^n\rangle\f$ */
    explicit AffineState(int n);

    /*! \brief Copy constructor*/
    AffineState(const AffineState& other);

    /*! \brief Controlled-\f$ Z \f$ gate*/
    void CZ(int j, int k);

    /*! \brief Controlled-\f$ X \f$ gate
     * \param j Control qubit
     * \param k Target qubit
     */
    void CX(int j, int k);

    /*! \brief Swap gate*/
    void SWAP(int j, int k);

    /*! \brief The gate \f$S = |0\rangle\langle0| + i|1\rangle\langle1|\f$ */
    void S(int j);

    /*! \brief The gate \f$S^\dagger = |0\rangle\langle0| -
     * i|1\rangle\langle1|\f$*/
    void SDG(int j);

    /*! \brief Hadamard gate */
    void H(int j);

    /*! \brief Pauli \f$ X \f$ gate*/
    void X(int j);

    /*! \brief Pauli \f$ Y \f$ gate*/
    void Y(int j);

    /*! \brief Pauli \f$ Z \f$ gate*/
    void Z(int j);

    /*! \brief Measure a single qubit
     * \param j Index of qubit to be measured
     * \param postselect Optional flag indicating whether postselection should
     * be used \param postselected_outcome Desired measurement outcome if
     * postselection is used. Throws error if postselection is chosen but is
     * impossible. \return Measurement outcome (0 or 1)
     */
    int MeasureZ(int j, bool postselect = false, int postselected_outcome = 0);

    /*! \brief Reset qubit to \f$|0\rangle\f$*/
    void Reset(int j);

    /*! \brief Measure all qubits
     * \return Vector of measurement outcomes
     */
    std::vector<int> MeasureAll() const;

    /*! \brief Sample the state repeatedly
     * \param nreps Number of samples
     * \return Map whose keys are vectors of measurement outcomes and whose
     * values are number of observations with each result
     */
    std::map<std::vector<int>, int> Sample(int nreps) const;

#ifdef USE_QPP
    /*!\brief Construct full state vector*/
    qpp::ket to_ket() const;
#endif // USE_QPP

    // Get member variables:

    /*! \brief Get number of qubits*/
    int n() const;

    /*! \brief Get value \f$ p \f$ such that global phase is \f$ e^{ip\pi/4}\f$
     */
    int phase() const;

    /*! \brief Get square matrix \f$ Q \f$ containing phases */
    mat_u_t Q() const;

    /*! \brief Get the generating matrix \f$ A \f$ of the affine space*/
    mat_u_t A() const;

    /*! \brief Get shift vector \f$ b \f$*/
    vec_u_t b() const;

    /*! \brief Get the so-called "principal index map." Entry \f$ i \f$
     * corresponds to the value \f$ p(i) \f$.*/
    std::vector<int> pivots() const;

    /*! \brief Get rank \f$ r \f$ of generating matrix \f$ A \f$*/
    int r() const;

    /*! \brief Print state in human-readable format*/
    friend std::ostream& operator<<(std::ostream& out, AffineState const& psi);

  private:
    int n_, phase_, r_;       // number of qubits, phase, and rank of A
    mat_u_t Q_, A_;           // Phase matrix and affine space generating matrix
    vec_u_t b_;               // Affine space shift vector
    std::vector<int> pivots_; // "Principal index map"

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

template <typename Derived>
using expr_t =
    typename std::decay<decltype(std::declval<Derived>().eval())>::type;

// Small helper functions
auto mod_2 = []() {
    // returns a unary lambda
    return [](unsigned x) { return x & 1; };
};

auto inline m2 = mod_2();

template <typename Derived>
[[nodiscard]] expr_t<Derived> ReduceMod2(const Eigen::MatrixBase<Derived>& A) {
    return A.unaryExpr(m2);
}

} /* namespace stab */

#endif /* STAB_AFFINE_STATE_H_ */
