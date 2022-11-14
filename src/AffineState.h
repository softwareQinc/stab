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

        void CZ(int a, int b);
        void CX(int a, int b);
        void SWAP(int a, int b);
        void S(int j);
        void SDG(int j);
        void H(int a);
        void X(int j);
        void Y(int j);
        void Z(int j);
        int MeasureZ(int j, bool postselect = false,
                     int postselected_outcome = 0);
        void Reset(int j);
        std::vector<int> MeasureAll() const;
        std::map<std::vector<int>, int> Sample(int nreps) const;
        Eigen::VectorXcd to_vec() const;

        // Get member variables:
        int n() const;

        int phase() const;

        mat_u_t Q() const;

        mat_u_t A() const;

        vec_u_t b() const;

        std::vector<int> pivots() const;

        int r() const;

        friend std::ostream &operator<<(std::ostream &out, AffineState const &psi);

    private:
        // State is represented by exp(i*pi*phase_/8)/sqrt(2^r_) \sum_{x \in
        // \{0,1\}^r_} i^{x^T Q_ x} \ket{Ax + b_ mod 2}, where A_ is n_\times r_ and
        // has rank r_ <= n_.
        int n_, phase_, r_;     // number of qubits
        mat_u_t Q_, A_;
        vec_u_t b_;
        std::vector<int> pivots_; // pivots_[i] == the row of A_ that has a pivot in column i

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
