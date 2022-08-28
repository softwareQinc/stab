#ifndef AFFINE_STATE_H_
#define AFFINE_STATE_H_

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <random>

namespace stab {
    class AffineState {
    public:

        explicit AffineState(int n);

        // Gates and measurements
        void CZ(int a, int b);

        void CX(int a, int b);

        void S(int j);

        void H(int a);

        void X(int j);

        void Y(int j);

        void Z(int j);

        int MeasureZ(int j); // Returns outcome and updates state
        void Reset(int j); // Resets qubit j to |0>

        void print_amplitudes();

        // Option to print state
        friend std::ostream &operator<<(std::ostream &out, AffineState const &psi);

    private:
        // State is represented by exp(i*pi*phase_/8)/sqrt(2^r_) \sum_{x \in
        // \{0,1\}^r_} i^{x^T Q_ x} \ket{Ax + b_ mod 2}, where A_ is n_\times r_ and
        // has rank r_ <= n_.
        int n_;     // number of qubits
        int phase_; // global phase_ is exp(i*pi*phase_/8)
        Eigen::MatrixXi Q_;
        Eigen::MatrixXi A_;
        Eigen::VectorXi b_;
        std::map<int, int>
                pivots_; // AKA "principal index map." Keys are columns, values are the
        // rows that contain pivots_ in those columns. Note that we are
        // using zero-indexing, so the smallest key (assuming the map
        // is nonempty) will always be 0, and the corresponding value
        // is the index of the row that has a pivot in column 0.
        // TODO: pivots can probably be changed to an std::unordered_map
        int r_; // Technically unnecessary since this is usually A_.cols() or rank of A_, but it is
        // handy to not have to declare it each time

        // Subroutines
        void FixFinalBit(int z);

        void ReduceGramRowCol(int c);

        void ReindexSubtColumn(int k, int c);

        void ReindexSwapColumns(int k, int c);

        void MakePrincipal(int c, int j);

        void ReselectPrincipalRow(int j, int c);

        void ZeroColumnElim(int c);

        void
        ReduceQ(); // Reduces mod 4 on the diagonal and mod 2 on the off-diagonal

        void print() const; // Only used for printing state at intermediate stages
        // of the calculation when diagnosing errors
    }; // class AffineState

// Small helper functions
    [[nodiscard]] Eigen::MatrixXi ReduceMod(const Eigen::MatrixXi &A, int modulus);
} // namespace stab

#endif // AFFINE_STATE_H_
