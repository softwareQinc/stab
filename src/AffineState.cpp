#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <vector>
#include <bitset>

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
    void AffineState::ReduceGramRowCol(int c) {
        int new_qcc = (4 + Q_(c, c) % 4) % 4; // Need to store this since it gets reduced mod 2 below
        ReduceMod(Q_.row(c), 2);
        ReduceMod(Q_.col(c), 2);
        Q_(c, c) = new_qcc;
    }

    void AffineState::ReindexSubtColumn(int k, int c) {
        // Applies the update Column k of A <--- Column k - Column c
        assert(c != k);
        A_.col(k) += A_.col(c);
        ReduceMod(A_.col(k), 2);
            
        int qcc = Q_(c, c);
        Q_.col(k) -= Q_.col(c);
        Q_.row(k) -= Q_.row(c);
        ReduceGramRowCol(k);
    }

    void AffineState::ReindexSwapColumns(int k, int c) {
        if (c != k) {
            A_.col(k).swap(A_.col(c));
            Q_.col(k).swap(Q_.col(c));
            Q_.row(k).swap(Q_.row(c));
            std::swap(pivots_.at(k), pivots_.at(c));
        }
    }

    void AffineState::MakePrincipal(int c, int j) {
        assert(A_(j, c) != 0);
        for (int k = 0; k < r_; k++) {
            if (k != c && A_(j, k) != 0) {
                ReindexSubtColumn(k, c);
            }
        }
        pivots_[c] = j;
    }

    void AffineState::ReselectPrincipalRow(int j, int c) {
        int j_star = -1; // Equivalent to j_star = 0 in the paper (but we need to
        // use -1 because of 0-indexing)
        for (int jj = 0; jj < A_.cols();
             ++jj) { //TODO: Choose optimal j_star (see paper)
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
        assert(r_ > 0);

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
        if (z == 1) { // ReduceMod is relatively expensive, so it makes sense to check whether z == 1
            Q_ += 2 * q.asDiagonal();
            ReduceMod(Q_.diagonal(), 4);
            b_ += a;
            ReduceMod(b_, 2);
            phase_ = (phase_ + 2 * u) % 8;
        }

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
        Q_(r_, r_) = 2 * b_(j);

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
            // zero columns, so this block would do nothing anyway
            Eigen::VectorXi a_j = A_.row(j).transpose();
            Eigen::VectorXi a_k = A_.row(k).transpose();
            Q_ += a_j * a_k.transpose() + a_k * a_j.transpose() +
                  2 * Eigen::MatrixXi(
                          (b_(k) * a_j.asDiagonal() + b_(j) * a_k.asDiagonal()));
            ReduceQ();
        }
        phase_ = (phase_ + 4 * b_(j) * b_(k)) % 8;
    }

    void AffineState::CX(int h, int j) {
        for (int i = 0; i < A_.cols(); ++i) {
            if (A_.col(i).isZero()) {
               std::cout << "Zero column entering CX subroutine\n";
               print();
            }
        }

        H(j);
        CZ(h, j);
        H(j);

        //// Step 1:
        //int c = -1;
        //for (int cc = 0; cc < r_; ++cc) {
        //    if (pivots_[cc] == j) {
        //        c = cc;
        //        break;
        //    }
        //}

        //// Step 2:
        //A_.row(j) += A_.row(h);
        //ReduceMod(A_.row(j), 2);

        //// Step 3:
        //b_(j) = (b_(j) + b_(h)) % 2;

        //// Step 4:
        //if (c != -1) {
        //    ReselectPrincipalRow(-1, c);
        //}

        for (int i = 0; i < A_.cols(); ++i) {
            if (A_.col(i).isZero()) {
                std::cout << "Zero column exiting CX on " << h << ", " << j
                          << "\n";
                print();
            }
        }
    }

    void AffineState::SWAP(int j, int k) {
        A_.row(j).swap(A_.row(k));
        b_.row(j).swap(b_.row(k));
    }

    void AffineState::S(int j) {
        bool is_pivot = false;
        int pivcol;
        for (int i = 0; i < r_; ++i) {
            if (pivots_[i] == j) {
                is_pivot = true;
                pivcol = i;
            }
        }

        if (is_pivot) {
            Q_(pivcol, pivcol) = (Q_(pivcol, pivcol) + (1 - 2 * b_(j)) + 4) % 4;
        } else {
            Eigen::VectorXi a_j = A_.row(j).transpose();
            Q_ += (1 - 2 * b_(j)) * a_j * a_j.transpose();
            ReduceQ();
        }
        phase_ = (phase_ + 2 * b_(j)) % 8;
    }

    void AffineState::X(int j) { b_(j) = (b_(j) + 1) % 2; }

    void AffineState::Z(int j) {
        phase_ = (phase_ + 4 * (b_(j))) % 8;
        Q_.diagonal() += 2 * A_.row(j);
        ReduceMod(Q_.diagonal(), 4);
    }

    void AffineState::Y(int j) {
        phase_ = (phase_ + 2) % 8;
        Z(j);
        X(j);
    }

    int AffineState::MeasureZ(int j) {
        // Check whether deterministic:
        if (A_.row(j).isZero()) {
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

    void AffineState::ReduceQ() {
        // Reduces Q mod 4 on the diagonal and mod 2 elsewhere
        Eigen::VectorXi qdiag = Q_.diagonal();
        ReduceMod(Q_, 2);
        ReduceMod(qdiag, 4);
        Q_.diagonal() = qdiag;
    }

    std::ostream &operator<<(std::ostream &out, AffineState const &psi) {
        out << "STATE IS GIVEN BY: \n";
        //out << "n = " << psi.n_ << '\n';
        out << "phase = exp(" << psi.phase_ << "*i*pi/8)\n";
        if (psi.r_ > 0) { // If psi is a computational basis state, no need to print this stuff
            out << "Q = " << std::endl << psi.Q_ << '\n';
            out << "A = " << std::endl << psi.A_ << '\n';
        }
        out << "b^T = " << psi.b_.transpose() << '\n';
        return out;
    }

    void AffineState::print() const { std::cout << '\n' << *this << std::endl; }

    void AffineState::print_amplitudes() {
        // Very naive way of printing the state, but this will mainly be used for debugging. We can always optimize it later if it's important.
        if (r_ > 15) {
            std::cout << "Too many amplitudes to print\n";
        } else {
            std::cout << "Global phase exp(" << phase_ << "*i*pi/8)\n";
            Eigen::VectorXi x;
            x.setZero(r_);
            for (int i = 0; i < pow(2, r_); ++i) {
                std::bitset<16> bs(i); // Get binary string
                for (int j = 0; j < r_; ++j) { // Cast binary string to x
                    x(j) = int(bs[j]);
                }

                Eigen::VectorXi ket = A_ * x + b_;
                ReduceMod(ket, 2);
                std::string rel_phase;
                switch ((x.transpose() * Q_ * x) % 4) {
                    case 0:
                        rel_phase = "  ";
                        break;
                    case 1:
                        rel_phase = " i";
                        break;
                    case 2:
                        rel_phase = " -";
                        break;
                    case 3:
                        rel_phase = "-i";
                        break;
                    default:
                        rel_phase = "abc";
                }

                std::cout << rel_phase << "|" << ket.transpose() << ">\n";
            }
        }
    }

} // namespace stab
