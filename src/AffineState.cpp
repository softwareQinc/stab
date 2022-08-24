#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"

namespace stab {
// Initialization:
    AffineState::AffineState(int n) : n_{n}, phase_{0}, r_{0} {
        Q_.setZero(0, 0);
        A_.setZero(n, 0);
        b_.setZero(n);
    }

// Subroutines:
    void AffineState::ReduceGramRowCol(int c) { // TODO: Optimize
        Q_(c, c) = (4 + Q_(c, c) % 4) % 4;
        for (int i = 0; i < Q_.cols(); i++) {
            if (i != c) {
                Q_(c, i) = (2 + Q_(c, i) % 2) % 2;
                Q_(i, c) = (2 + Q_(i, c) % 2) % 2;
            }
        }
    }

    void AffineState::ReindexSubtColumn(int k, int c) {
        // Applies the update Column k of A <--- Column k - Column c
        if (c != k) {
            for (int j = 0; j < n_; j++) {
                A_(j, k) = (A_(j, k) + A_(j, c)) % 2;
            }
            int qcc = Q_(c, c);
            for (int h = 0; h < Q_.cols(); h++) {
                if (Q_(h, c) != 0) {
                    //Q_(h, k) -= Q_(h, c);
                    Q_(h, k) += Q_(h, c); // TODO: Check whether +/- is correct
                }
                if (Q_(c, h) != 0) {
                    //Q_(k, h) -= Q_(c, h);
                    Q_(k, h) += Q_(c, h); // TODO: Check whether +/- is correct
                }
            }
            Q_(k, k) += qcc; // TODO: Check whether correct
            ReduceGramRowCol(k);
            ReduceQ();
        }
    }

    void AffineState::ReindexSwapColumns(int k, int c) {
        if (c != k) {
            int swaptmp;
            for (int j = 0; j < n_; j++) {
                if (A_(j, c) != 0 || A_(j, k) != 0) { // TODO: Optimize swap
                    swaptmp = A_(j, c);
                    A_(j, c) = A_(j, k);
                    A_(j, k) = swaptmp;
                }
            }
            for (int j = 0; j < Q_.cols(); j++) {
                if (Q_(j, c) != 0 || Q_(j, k) != 0) { // TODO: Optimize swap
                    swaptmp = Q_(j, c);
                    Q_(j, c) = Q_(j, k);
                    Q_(j, k) = swaptmp;
                }
            }
            for (int j = 0; j < Q_.cols(); j++) {
                if (Q_(c, j) != 0 || Q_(k, j) != 0) { // TODO: Optimize swap
                    swaptmp = Q_(c, j);
                    Q_(c, j) = Q_(k, j);
                    Q_(k, j) = swaptmp;
                }
            }
            ReduceQ();
            std::swap(pivots_.at(k), pivots_.at(c));
        }
    }

    void AffineState::MakePrincipal(int c, int j) {
        if (A_(j, c) == 0) {
            std::cerr << "Error......." << std::endl;
            return;
        } else {
            for (int k = 0; k < r_; k++) {
                if (k != c && A_(j, k) != 0) {
                    ReindexSubtColumn(k, c);
                }
            }
            pivots_[c] = j;
        }
    }

    void AffineState::ReselectPrincipalRow(int j, int c) {
        int j_star = -1; // Equivalent to j_star = 0 in the paper (but we need to
        // use -1 because of 0-indexing)
        for (int jj = 0; jj < A_.cols();
             jj++) { // TODO: Rewrite to find optimal j_star (see paper)
            if (A_(jj, c) != 0 && jj != j) {
                j_star = jj;
                break;
            }
        }
        if (j_star != -1) {
            MakePrincipal(c, j_star);
        } // (else, do nothing since reselecting row is impossible)
        return;
    }

    void AffineState::FixFinalBit(int z) {
        // TODO: Confirm that whenever this function is called, A_ will have r_ columns and rank = r_

        // Step 1:
        Eigen::VectorXi a = A_(Eigen::all, r_ - 1); // r_ - 1 because of 0-indexing
        Eigen::VectorXi q;
        if (r_ > 1) { // If r == 1 then the line below throws an error due to r_ - 2
            q = Q_(Eigen::seq(0, r_ - 2), r_ - 1); 
        }
        int u = Q_(r_ - 1, r_ - 1) % 4;
        // Step 2:
        A_.conservativeResize(Eigen::NoChange, r_ - 1); // Keep only first r_ - 1 columns
        Q_.conservativeResize(r_ - 1, r_ - 1);

        // Step 3:
        Q_ += 2 * z * q.asDiagonal();
        ReduceQ(); // TODO: Possibly unnecessary
        b_ += z * a;
        ReduceVectorMod(b_, 2);
        phase_ = (phase_ + 2 * z * u) % 8;

        // Step 4:
        pivots_.erase(r_ - 1);
        --r_;

        return;
    }

    void AffineState::ZeroColumnElim(int c) {
        // Recall: Matrix A_ has rank r_ and r_+1 columns at this stage

        // Step 1:
        ReindexSwapColumns(c, r_); // Use r_ instead of r_+1 here because of 0-indexing

        // Step 2:
        Eigen::VectorXi q = Q_(Eigen::seq(0, r_ - 1), r_); // Step 2.
        int u = Q_(r_, r_) % 4;

        // Step 3:
        A_.conservativeResize(Eigen::NoChange, r_); // A^(2) from paper
        Q_.conservativeResize(r_, r_); // Q^(2) from paper

        // Step 3.5:
        pivots_.erase(r_); // <-- Both parts of the if statement in Step 4 require us to do this

        // Step 4
        if (u % 2 == 1) {
            Q_ += (2 - u) * q * q.transpose(); // TODO: Confirm that the "u-2" from the paper is an error
            ReduceQ();
            phase_ = (phase_ - u + 2) % 8;
            return;
        } else {
            // TODO: Check whether q == \vec{0}? For now assume not
            int ell;
            for (int i = 0; i < r_; i++) {
                if (q(i) != 0) {
                    ell = i; // Because q is guaranteed to be nonzero, we are guaranteed to initialize ell eventually
                    break;
                }
            }
            for (int k = 0; k < r_; k++) {
                if (q(k) != 0 && k != ell) {
                    ReindexSubtColumn(k, ell);
                }
            }
            ReindexSwapColumns(r_ - 1, ell); // r_-1 because of 0-indexing. This step also produces A^(3) and Q^(3)
            // At this point, A_ and Q_ should have r_ columns each. Also, since we removed the zero column in Step 3, the rank of A_ should be r_
            FixFinalBit(u / 2);
            return;
        }
    }

// GATES:
    void AffineState::H(int j) {
        ReduceQ(); // TODO: Removing this probably is fine

        // Step 1: Find c.
        int c = -1; // c = -1 will indicate that row j is not a pivot
        for (int column_i = 0; column_i < r_; column_i++) {
            if (pivots_[column_i] == j) {
                c = column_i;
                break;
            }
        }
        
        // Step 2. (Using -1 instead of 0 because of 0-indexing.)
        if (c > -1) {
            ReselectPrincipalRow(j, c);
            if (pivots_[c] != j) { // Check if reselection was successful
                c = -1; // Indicates that row j is not/no longer a pivot
            }
        }

        // Step 3:
        Eigen::VectorXi atilde;
        atilde.setZero(r_ + 1);
        atilde(Eigen::seq(0, r_ - 1)) =
                A_(j, Eigen::all)
                        .transpose(); // Set atilde = [row j of A_ | 0 ]^T \in \{0,1\}^{r+1}

        // Step 4:
        A_.row(j).setZero();
        A_.conservativeResize(
                Eigen::NoChange,
                r_ + 1); // Add column, and (next two lines) set new column to e_j
        A_.col(r_).setZero();
        A_(j, r_) = 1;
        pivots_[r_] = j;

        // Step 5:
        Q_.conservativeResize(r_ + 1, r_ + 1);
        Q_.row(r_) = atilde.transpose();
        Q_.col(r_) = atilde.transpose();
        Q_(r_, r_) = 2 * b_(j); // Already reduced mod 4
        ReduceQ(); // TODO: Probably unnecessary

        // Step 6:
        b_(j) = 0;
        // Step 7:
        if (c > -1) { // Case 2 from Alex's notes
            ZeroColumnElim(c);
        } else { // Case 1 from Alex's notes
            ++r_;
        }
    }

    void AffineState::CZ(int j, int k) {
        if (A_.cols() > 0) { // If psi is a computational basis state then A_ has
            // zero columns, so the next three lines would do nothing anyway
            Eigen::VectorXi a_j = A_(j, Eigen::all).transpose();
            Eigen::VectorXi a_k = A_(k, Eigen::all).transpose();
            Q_ += a_j * a_k.transpose() + a_k * a_j.transpose() +
                  2 * Eigen::MatrixXi(
                          (b_(k) * a_j.asDiagonal() + b_(j) * a_k.asDiagonal()));
            ReduceQ();
        }
        phase_ = (phase_ + 4 * b_(j) * b_(k)) % 8;
    }

    void AffineState::CX(int j, int k) { // TODO: Use non-naive method.
        H(k);
        CZ(j, k);
        H(k);
    }

    void AffineState::S(int j) {
        Eigen::VectorXi a_j = A_(j, Eigen::all).transpose();
        Q_ += (1 - 2 * b_(j)) * a_j * a_j.transpose();
        ReduceQ();
        phase_ = (phase_ + 2 * b_(j)) % 8;
    }

    void AffineState::X(int j) { b_(j) = (b_(j) + 1) % 2; }

    void AffineState::Z(int j) {
        phase_ = (phase_ + 4 * (b_(j))) % 8;
        Q_ += 2 * Eigen::MatrixXi(A_(j, Eigen::all).asDiagonal());
        ReduceQ();
    }

    void AffineState::Y(int j) {
        phase_ = (phase_ + 2) % 8;
        Z(j);
        X(j);
    }

    int AffineState::MeasureZ(int j) {
        // Check whether deterministic:
        if (A_(j, Eigen::all).isZero()) {
            return b_(j);
        } else {
            int beta = random_bit();
            int k;
            for (int i = 0; i < r_; i++) { // TODO: Optimize
                if (A_(j, i) != 0) {
                    k = i;
                    break;
                }
            }
            ReindexSwapColumns(k, r_ - 1);
            MakePrincipal(r_ - 1, j);
            FixFinalBit(beta);
            return beta;
        }
    }

    void AffineState::Reset(int j) {
        int tmp = MeasureZ(j);
        if (tmp == 1) {
            X(j);
        }
    }

    void ReduceMatrixMod(Eigen::MatrixXi &M,
                         int modulus) { // TODO: This is a really silly way of doing
        // it. Should find something better
        M = M.array() - (modulus * (M.array() / modulus));
        M = M.array() + modulus;
        M = M.array() - (modulus * (M.array() / modulus));
    }

    void ReduceVectorMod(Eigen::VectorXi &v,
                         int modulus) { // TODO: This is a really silly way of doing
        // it. Should find something better
        v = v.array() - (modulus * (v.array() / modulus));
        v = v.array() + modulus;
        v = v.array() - (modulus * (v.array() / modulus));
    }

    void AffineState::ReduceQ() { // TODO: Optimize
        // Need to reduce off-diagonal elements modulo 2 and diagonal elements
        // modulo 4
        /*std::cout << "Reducing Q..." << std::endl << std::endl;
        std::cout << "Q = " << std::endl << Q_ << std::endl;*/
        Eigen::MatrixXi Qdiag = Q_.diagonal().asDiagonal();
        Q_ = Q_ - Qdiag;
        ReduceMatrixMod(Q_, 2);
        ReduceMatrixMod(Qdiag, 4);
        Q_ += Qdiag;
        /*std::cout << "Now, Q = " << std::endl << Q_ << std::endl;*/
    }

    std::ostream &operator<<(std::ostream &out, AffineState const &psi) {
        out << "STATE IS GIVEN BY: \n";
        //out << "n = " << psi.n_ << '\n';
        //out << "r = " << psi.r_ << '\n';
        out << "phase = " << psi.phase_ << '\n';
        out << "Q = " << std::endl << psi.Q_ << '\n';
        out << "A = " << std::endl << psi.A_ << '\n';
        out << "b^T = " << psi.b_.transpose() << '\n';
        /*std::cout << "pivots = ";
        for (const auto &p: psi.pivots_) {
            out << "(" << p.first << ", " << p.second << "), ";
        }*/
        return out;
    }

    void AffineState::print() const { std::cout << '\n' << *this << std::endl; }

    void AffineState::print_amplitudes() {
        // TODO: It would be nice to have some way to print out each ket and amplitude
    }

} // namespace stab
