/**
\file AffineState.h
\brief All functions needed for affine state-based simulation described in
arXiv:2109.08629
\author softwareQ
*/

#ifndef STAB_AFFINE_STATE_H_
#define STAB_AFFINE_STATE_H_

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

///< maximum number of qubits allowed in qpp::AffineState::to_ket()
#define MAX_QUBITS_STATE_VECTOR 16

namespace stab {
class AffineState {
  public:
    /*! \brief Create the state \f$|0^n\rangle\f$ */
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

    /*! \brief Resets qubit to \f$|0\rangle\f$*/
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

    /*!\brief Construct full statevector*/
    Eigen::VectorXcd to_ket() const;

    // Get member variables:

    /*! \brief Get number of qubits*/
    int n() const;

    /*! \brief Get value \f$ p \f$ such that global phase is \f$ e^{ip\pi/4}\f$
     */
    int phase() const;

    /*! \brief Get \f$ Q \f$ matrix containing phases */
    mat_u_t Q() const;

    /*! \brief Get \f$ A \f$ generating the affine space*/
    mat_u_t A() const;

    /*! \brief Get offset vector \f$ b \f$*/
    vec_u_t b() const;

    /*! \brief Get the principal index map. Entry \f$ i \f$ corresponds to the
     * value \f$ p(i) \f$.*/
    std::vector<int> pivots() const;

    /*! \brief Get rank \f$ r \f$ of generating matrix*/
    int r() const;

    /*! \brief Print state in human-readable format*/
    friend std::ostream& operator<<(std::ostream& out, AffineState const& psi);

  private:
    int n_, phase_, r_; // number of qubits
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
