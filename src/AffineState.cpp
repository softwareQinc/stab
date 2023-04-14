#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <bitset>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"

namespace stab {
// constructor (initialization):
    AffineState::AffineState(int n) : n_{n}, phase_{0}, r_{0} {
        pivots_ = std::vector<int>(); 
        pivots_.reserve(n + 1);
        Q_.setZero(n + 1, n + 1);
        A_.setZero(n, n + 1);
        b_.setZero(n);
    }

    // copy
    AffineState::AffineState(const AffineState &other) {
        n_ = other.n_;
        phase_ = other.phase_;
        r_ = other.r_;
        Q_ = other.Q_;
        A_ = other.A_;
        b_ = other.b_;
        pivots_ = other.pivots_;
    }

    int AffineState::n() const { return n_; }

    int AffineState::phase() const { return phase_; }

    mat_u_t AffineState::Q() const { return Q_.block(0, 0, r_, r_); }

    mat_u_t AffineState::A() const { return A_.block(0, 0, n_, r_); }

    vec_u_t AffineState::b() const { return b_; }

    std::vector<int> AffineState::pivots() const { return pivots_; }

    int AffineState::r() const { return r_; }

// Subroutines:
    std::vector<int> AffineState::A_col_nonzeros(int col) {
        // Returns location of ones in column "col" of A
        std::vector<int> ones;
        ones.reserve(n_);
        for (int row = 0; row < n_; ++row) {
            if (A_(row, col) == 1)
                ones.push_back(row);
        }
        return ones;
    }

    std::vector<int> AffineState::Q_nonzeros() {
        // Returns location of nonzeros in last col of Q, ignoring the diagonal entry
        std::vector<int> ones;
        ones.reserve(r_ - 1);
        for (int row = 0; row < r_ - 1; ++row) {
            if (Q_(row, r_ - 1) != 0)
                ones.push_back(row);
        }
        return ones;
    }

    std::vector<int> AffineState::A_row_nonzeros(int row) {
        std::vector<int> ones;
        ones.reserve(r_);
        for (int col = 0; col < r_; ++col) {
            if (A_(row, col) != 0)
                ones.push_back(col);
        }
        return ones;
    }

    void AffineState::ReduceGramRowCol(int c) {
        int new_qcc = Q_(c, c) % 4; // Need to store since gets reduced mod 2 below
        Q_.row(c) = ReduceMod2(Q_.row(c));
        Q_.col(c) = ReduceMod2(Q_.col(c));
        Q_(c, c) = new_qcc;
    }

    void AffineState::ReindexSubtColumn(int k, int c, std::vector<int> col_c_nonzeros) {
        // Applies the update Column k of A <--- Column k - Column c
        assert(c != k);

        for (int row : col_c_nonzeros) A_(row, k) ^= 1;

        // TODO: Optimize Q part? Worst-case complexity won't change
        Q_.col(k) += Q_.col(c);
        Q_.row(k) += Q_.row(c);
        ReduceGramRowCol(k);
    }

    void AffineState::ReindexSwapColumns(int k) {
        if (k != r_ - 1) {
            A_.col(k).swap(A_.col(r_ - 1));
            Q_.col(k).swap(Q_.col(r_ - 1));
            Q_.row(k).swap(Q_.row(r_ - 1));
            std::swap(pivots_[k], pivots_[r_ - 1]);
        }
    }

    void AffineState::MakePrincipal(int c, int j) {
        assert(A_(j, c) != 0);
        std::vector<int> col_c_nonzeros = A_col_nonzeros(c);
        for (int k = 0; k < r_; ++k) {
            if (A_(j, k) != 0 && k != c) {
                ReindexSubtColumn(k, c, col_c_nonzeros);
            }
        }
        pivots_[c] = j;
    }

    bool AffineState::ReselectPrincipalRow(int j, int c) {
        // Performs update if possible and returns flag indicating success or failure
        for (int j_star = 0; j_star < n_; ++j_star) {
            if (A_(j_star, c) != 0 && j_star != j) {
                MakePrincipal(c, j_star);
                return true;
            }
        }
        return false;
    }

    void AffineState::FixFinalBit(int z) {
        assert(r_ != 0);

        // Step 1:
        std::vector<int> a_ones = A_col_nonzeros(r_ - 1);
        std::vector<int> q_ones = Q_nonzeros();
        int u = Q_(r_ - 1, r_ - 1) % 4;

        // Step 2:
        for (int row : a_ones) A_(row, r_ - 1) = 0; 
        
        for (int rowcol : q_ones) {
            Q_(rowcol, r_ - 1) = 0;
            Q_(r_ - 1, rowcol) = 0;
        }
        Q_(r_ - 1, r_ - 1) = 0;

        // Step 3:
        if (z == 1) {
            for (int i : q_ones) Q_(i, i) = (Q_(i, i) + 2) % 4;
            for (int row : a_ones) b_(row) ^= 1;
            phase_ = (phase_ + 2 * u) % 8;
        }

        // Step 4:
        pivots_.pop_back();
        --r_;
    }

    void AffineState::ZeroColumnElim(int c) {
        // Step 1:
        ReindexSwapColumns(c);

        // Step 2:
        std::vector<int> q_ones = Q_nonzeros();
        int u = Q_(r_ - 1, r_ - 1) % 4;

        // Step 3:
        for (int rowcol : q_ones) {
            Q_(rowcol, r_-1) = 0;
            Q_(r_-1, rowcol) = 0;
        }
        Q_(r_-1, r_-1) = 0;

        // Step 3.5:
        pivots_.pop_back();
        --r_;

        // Step 4
        if (u % 2 == 1) {
            for (int row : q_ones) {
                for (int col : q_ones) {
                    if (row == col) {
                        Q_(row, col) = (Q_(row, col) + 2 + u) % 4;
                    } else {
                        Q_(row, col) ^= 1;
                    }
                }
            }
            phase_ = (phase_ - u + 2) % 8;
        } else { // u is even
            assert(q_ones.size() != 0);

            int ell = q_ones[0];
            std::vector<int> col_ell_nonzeros = A_col_nonzeros(ell);

            for (int col : q_ones) {
                if (col != ell) ReindexSubtColumn(col, ell, col_ell_nonzeros); // TODO: Check if ok w/r/t r_
            }

            ReindexSwapColumns(ell); // r_-1 because of 0-indexing. This step also produces A^(3) and Q^(3)
            // At this point, A_ and Q_ should have r_ columns each. Also, since we removed the zero column in Step 3, the rank of A_ should be r_
            FixFinalBit(u / 2);
        }
    }

    int AffineState::piv_col(int row_number) {
        // Given a row number, return the index of the column j for in which that row has a pivot, and return -1 otherwise
        // TODO: Optimize
        for (int i = 0; i < r_; ++i) {
            if (pivots_[i] == row_number) {
                return i;
            }
        }
        return -1;
    }

// GATES:
    void AffineState::H(int j) {
        // Step 1: Find c.
        int c = piv_col(j); // c = -1 indicates that row j is not a pivot
        // Step 2. (Using -1 instead of 0 because of 0-indexing.)
        if (c != -1) {
            if (ReselectPrincipalRow(j, c)) {
                c = -1;
            }
        }

        // Step 3:
        std::vector<int> atilde_ones = A_row_nonzeros(j);

        // Step 4:
        for (int col : atilde_ones) {
            A_(j, col) = 0;
        }
        
        A_(j, r_) = 1;
        pivots_.push_back(j);
        ++r_;

        // Step 5:
        for (int rowcol : atilde_ones) {
            Q_(rowcol, r_-1) = 1;
            Q_(r_-1, rowcol) = 1;
        }
        Q_(r_-1, r_-1) = 2 * b_(j);

        // Step 6:
        b_(j) = 0;

        // Step 7:
        if (c > -1) ZeroColumnElim(c);
    }

    void AffineState::CZ(int j, int k) {
        std::vector<int> Aj_ones = A_row_nonzeros(j);
        std::vector<int> Ak_ones = A_row_nonzeros(k);

        for (int jj : Aj_ones) {
            for (int kk : Ak_ones) {
                if (jj == kk) {
                    Q_(jj, kk) = (Q_(jj, kk) + 2) % 4;
                } else {
                    Q_(jj, kk) ^= 1;
                    Q_(kk, jj) ^= 1;
                }
            }
            Q_(jj, jj) = (Q_(jj, jj) + 2 * b_(k)) % 4;
        }

        for (int kk : Ak_ones) Q_(kk, kk) = (Q_(kk, kk) + 2 * b_(j)) % 4;

        phase_ = (phase_ + 4 * b_(j) * b_(k)) % 8;
    }

    void AffineState::CX(int h, int j) {
        int c = piv_col(j);
        for (int col : A_row_nonzeros(h)) A_(j, col) ^= 1;
        b_(j) ^= b_(h);
        if (c != -1) ReselectPrincipalRow(-1, c);
    }

    void AffineState::SWAP(int j, int k) {
        // TODO: Use better way to update pivots_
        A_.row(j).swap(A_.row(k));
        b_.row(j).swap(b_.row(k));

        auto it_j = std::find(pivots_.begin(), pivots_.end(), j);
        auto it_k = std::find(pivots_.begin(), pivots_.end(), k);
        if (it_j != pivots_.end()) *it_j = k;
        if (it_k != pivots_.end()) *it_k = j;
    }

    void AffineState::S_or_SDG(int j, bool dg) {
        // dg = true means we apply S^\dagger, dg = false means we apply S
        int sign = 1 - 2 * int(dg);  // = +1 for S and -1 for S^\dagger

        std::vector<int> Aj_ones = A_row_nonzeros(j);
        for (int row : Aj_ones) {
            for (int col : Aj_ones) {
                if (row == col) {
                    Q_(row, col) += 4 + sign * (1 - 2 * b_(j));
                    Q_(row, col) %= 4;
                } else {
                    Q_(row, col) ^= 1;
                }
            }
        }
        phase_ = (phase_ + sign * 2 * b_(j) + 8) % 8;
    }

    void AffineState::S(int j) { S_or_SDG(j, false); }

    void AffineState::SDG(int j) { S_or_SDG(j, true); }

    void AffineState::X(int j) { b_(j) ^= 1; }

    void AffineState::Z(int j) {
        phase_ = (phase_ + 4 * (b_(j))) % 8;
        for (int i = 0; i < r_; ++i) {
            if (A_(j, i) != 0) {
                Q_(i, i) = (Q_(i, i) + 2) % 4;
            }
        }
    }

    void AffineState::Y(int j) {
        phase_ = (phase_ + 2) % 8;
        Z(j);
        X(j);
    }

    int AffineState::MeasureZ(int j, bool postselect,
                              int postselected_outcome) {
        // Default parameters are postselect=false, postselected_outcome=0
        if (postselected_outcome != 0 && postselected_outcome != 1) {
            throw std::logic_error("Measurement outcome must be 0 or 1");
        }
        if (A_.row(j).isZero()) { // Deterministic case
            if (postselect && b_(j) != postselected_outcome) {
                throw std::logic_error("Postselection impossible");
            } else {
                return b_(j);
            }
        } else { // Random case
            int beta;
            if (postselect) {
                beta = postselected_outcome;
            } else {
                beta = random_bit();
            }
            for (int k = 0; k < r_; k++) {
                if (A_(j, k) != 0) {
                    ReindexSwapColumns(k);
                    break;
                }
            }
            MakePrincipal(r_ - 1, j);
            FixFinalBit(beta^b_(j));
            return beta;
        }
    }

    std::vector<int> AffineState::MeasureAll() const {
        vec_u_t x(r_);
        for (int i = 0; i < r_; ++i) x(i) = random_bit();

        vec_u_t tmp = ReduceMod2(A_.block(0,0,n_,r_) * x + b_);

        std::vector<int> result(n_);
        for (int i = 0; i < n_; ++i) result[i] = tmp(i);

        return result;
    }

    std::map<std::vector<int>, int>
    AffineState::Sample(int nreps) const {
        std::map<std::vector<int>, int> results;
        for (int repnumber = 0; repnumber < nreps; ++repnumber) {
            ++results[MeasureAll()];
        }
        return results;
    }

    void AffineState::Reset(int j) {
        int tmp = MeasureZ(j);
        if (tmp == 1) X(j);
    }

    Eigen::VectorXcd AffineState::to_vec() const {
        // SUPER ugly way to get the statevector. Only needed for unit tests.
        if (n_ > 16) {
            throw std::logic_error("Maximum number of qubits for statevector "
                                   "representation is 15");
        }
        std::complex<double> pi = std::atan(1.0) * 4.0;
        std::complex<double> i(0, 1);


        Eigen::VectorXcd vec; // This will be the statevector
        vec.setZero(int(pow(2, n_)));
        for (int k = 0; k < pow(2, r_); ++k) {
            // For each x \in \{0,1\}^ncols, we now add the corresponding term to vec
            // First, let x = vector whose entries are the binary digits of k:
            vec_u_t x(r_);
            std::bitset<16> bs(k);             // Get binary string
            for (int j = 0; j < r_; ++j) { // Cast binary string to x
                x(j) = int(bs[j]);
            }

            // Figure out which basis state x results in:
            vec_u_t ket = ReduceMod2(A_.block(0,0,n_,r_) * x + b_);
            // "ket" is the binary representation of some number basis_state_number. We need to add the correct amplitude into vec[basis_state_number]. Note that AffineState puts the zero-th qubit on the left, but most people put it on the right, so we also make this conversion.

            int basis_state_number = 0;
            for (int j = 0; j < n_; ++j) {
                basis_state_number += int(ket(j) * pow(2, n_ - 1 - j)); // Conversion
            }

            std::complex<double> phase = (2 * (x.transpose() * Q_.block(0,0,r_,r_) * x)[0] + phase_);
            vec[basis_state_number] += std::exp(i * pi * phase / 4.0);
        }
        return vec / pow(2, r_ / 2.0);
    }

    std::ostream &operator<<(std::ostream &out, AffineState const &psi) {
        out << "\nQUADRATIC FORM REPRESENTATION: \n";
        out << "n = " << psi.n_ << "\n";
        out << "phase = exp(" << psi.phase_ << "*i*pi/4)\n";
        out << "r = " << psi.r_ << "\n";
        if (psi.r_ > 0) { // If psi is a computational basis state, no need to print this stuff
            out << "Q = " << std::endl << psi.Q_.block(0,0,psi.r_,psi.r_) << '\n';
            out << "A = " << std::endl << psi.A_.block(0,0,psi.n_,psi.r_) << '\n';
        }
        out << "b^T = " << psi.b_.transpose() << '\n';
        out << "pivots = ";
        for (int i = 0; i < psi.r_; ++i) {
            out << "(" << i << ", " << psi.pivots_[i] << ") ";
        }
        return out;
    }

} // namespace stab
