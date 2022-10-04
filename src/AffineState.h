#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>

#include <qpp/qpp.h>

using block_t = Eigen::Block<Eigen::MatrixXi>;

namespace stab {
    class AffineState {
    public:

        explicit AffineState(int n);

        // Unitary gates
        void CZ(int a, int b);

        void CX(int a, int b);

        void SWAP(int a, int b);

        void S(int j);

        void SDG(int j);

        void H(int a);

        void X(int j);

        void Y(int j);

        void Z(int j);

        // Nonunitary operations
        int MeasureZ(int j, bool postselect = false,
                     int postselected_outcome = 0);
        void Reset(int j); // Resets qubit j to |0>

        std::vector<int> MeasureAll() const;

        std::map<std::vector<int>, int> Sample(int nreps) const;

        Eigen::VectorXcd to_vec() const;

        // Get member variables:
        int n() const;

        int phase() const;

        Eigen::MatrixXi Q() const;

        Eigen::MatrixXi A() const;

        Eigen::VectorXi b() const;

        std::unordered_map<int, int> pivots() const;

        int r() const;

        friend std::ostream &operator<<(std::ostream &out, AffineState const &psi);

    private:
        // State is represented by exp(i*pi*phase_/8)/sqrt(2^r_) \sum_{x \in
        // \{0,1\}^r_} i^{x^T Q_ x} \ket{Ax + b_ mod 2}, where A_ is n_\times r_ and
        // has rank r_ <= n_.
        int n_;     // number of qubits
        int phase_; // global phase_ is exp(i*pi*phase_/4)
        Eigen::MatrixXi Qmaster_;
        Eigen::MatrixXi Amaster_;
        std::unique_ptr<block_t> Q_; // Points to active block of Qmaster_ (top-left corner)
        std::unique_ptr<block_t> A_; // Points to active block of A_ (leftmost columns)
        Eigen::VectorXi b_;
        std::unordered_map<int, int>
                pivots_; // AKA "principal index map." Keys are columns, values are the
        // rows that contain pivots_ in those columns. Note that we are
        // using zero-indexing, so the smallest key (assuming the map
        // is nonempty) will always be 0, and the corresponding value
        // is the index of the row that has a pivot in column 0.
        // TODO: pivots can probably be changed to an std::unordered_map
        int r_; // Technically unnecessary since this is usually the rank of A_, but it is
        // handy to not have to declare it each time

        // Subroutines
        void FixFinalBit(int z);

        void ReduceGramRowCol(int c);

        void ReindexSubtColumn(int k, int c);

        void ReindexSwapColumns(int k, int c);

        void MakePrincipal(int c, int j);

        bool ReselectPrincipalRow(int j, int c);

        void ZeroColumnElim(int c);

        void
        ReduceQ(); // Reduces mod 4 on the diagonal and mod 2 on the off-diagonal

        void S_or_SDG(int j, bool dg);

        int piv_col(int row_number);

        void print() const; // Only used for printing state at intermediate stages
        // of the calculation when diagnosing errors
    }; // class AffineState

    template<typename Derived>
    using expr_t =
            typename std::decay<decltype(std::declval<Derived>().eval())>::type;

    // Small helper functions
    template<typename Derived>
    [[nodiscard]] expr_t<Derived> ReduceMod(const Eigen::MatrixBase<Derived> &A, int modulus) {
        // Reduces matrix modulo "modulus," and ensures entries are all nonnegative
        assert(modulus > 0);
        // define the mod_p as a lambda taking an int (modulo) and returning a
        // unary lambda that computes x -> x mod p \in {0,...,p-1};
        auto mod_p = [](int p) {
            // returns a unary lambda
            return [p](int x) { return ((x % p) + p) % p; };
        };
        return A.unaryExpr(mod_p(modulus));
    }

} // namespace stab
