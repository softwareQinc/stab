#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>

#include <qpp/qpp.h>

using mat_b_t = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
using vec_b_t = Eigen::Vector<bool, Eigen::Dynamic>;
using block_t = Eigen::Block<mat_b_t>;
using subvec_t = Eigen::VectorBlock<vec_b_t>;

namespace stab {
    class AffineState {
    public:

        explicit AffineState(int n);

        AffineState(const AffineState &other);

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
        bool MeasureZ(int j, bool postselect = false,
                     bool postselected_outcome = false);

        void Reset(int j); // Resets qubit j to |0>

        std::vector<int> MeasureAll() const;

        std::map<std::vector<int>, int> Sample(int nreps) const;

        Eigen::VectorXcd to_vec() const;

        // Get member variables:
        int n() const;

        int phase() const;

        mat_b_t Q() const;

        mat_b_t A() const;

        vec_b_t b() const;

        std::unordered_map<int, int> pivots() const;

        int r() const;

        friend std::ostream &operator<<(std::ostream &out, AffineState const &psi);

    private:
        // State is represented by exp(i*pi*phase_/8)/sqrt(2^r_) \sum_{x \in
        // \{0,1\}^r_} i^{x^T Q_ x} \ket{Ax + b_ mod 2}, where A_ is n_\times r_ and
        // has rank r_ <= n_.
        int n_;     // number of qubits
        int phase_; // global phase_ is exp(i*pi*phase_/4)

        mat_b_t Qmaster_;
        mat_b_t Amaster_;
        vec_b_t Qonesmaster_;
        vec_b_t Qtwosmaster_;
        std::unique_ptr<block_t> Q_; // Points to active block of Qmaster_ (top-left corner)
        std::unique_ptr<block_t> A_; // Points to active block of A_ (leftmost columns)
        std::unique_ptr<subvec_t> Qones_; // ...
        std::unique_ptr<subvec_t> Qtwos_; // ...
        vec_b_t b_;
        std::unordered_map<int, int>
                pivots_; // AKA "principal index map." Keys are columns, values are the
        // rows that contain pivots_ in those columns. Note that we are
        // using zero-indexing, so the smallest key (assuming the map
        // is nonempty) will always be 0, and the corresponding value
        // is the index of the row that has a pivot in column 0.
        int r_;

        // Subroutines
        void DiagAdd(int index, int addend);

        void FixFinalBit(int z);

        void ReindexSubtColumn(int k, int c);

        void ReindexSwapColumns(int k, int c);

        void MakePrincipal(int c, int j);

        bool ReselectPrincipalRow(int j, int c);

        void ZeroColumnElim(int c);

        void S_or_SDG(int j, bool dg);

        int piv_col(int row_number);

        void print() const; // Only used for printing state at intermediate stages
        // of the calculation when diagnosing errors
    }; // class AffineState

    template<typename Derived>
    using expr_t =
            typename std::decay<decltype(std::declval<Derived>().eval())>::type;

    // Small helper function
    // Eigen uses "+" as the "OR" operator, but we need XOR, so we will use this instead
    auto bool_xor = [](bool x, bool y) { return x != y; };

} // namespace stab
