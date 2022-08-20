#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "AffineState.h"
#include "random.h"

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
    if (c != k) {
        for (int j = 0; j < n_; j++) {
            if (A_(j, c) != 0) {
                A_(j, k) = (A_(j, k) + 1) % 2;
            }
        }
        for (int h = 0; h < Q_.cols(); h++) {
            if (Q_(h, c) != 0) {
                Q_(h, k) -= Q_(h, c);
            }
            if (Q_(c, h) != 0) {
                Q_(k, h) -= Q_(c, h);
            }
        }
        ReduceGramRowCol(k);
        // ReduceQ();
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
    for (int jj = 0; jj < r_;
         jj++) { // TODO: Rewrite to find optimal j_star (see paper)
        if (A_(jj, c) != 0) {
            j_star = jj;
            break;
        }
    }
    if (j_star != -1) {
        MakePrincipal(c, j_star);
        return;
    } else { // TODO: Don't think is an error is necessary here, but should
             // doublecheck eventually.
        return;
    }
}

void AffineState::FixFinalBit(int z) {
    std::cout << "Entering FixFinalBit";
    print();
    int rtemp =
        A_.cols(); // Define this because sometimes this function gets called
                   // from the ZeroColumnElim function, in which r_ might
                   // temporarily not be out of sync with A_'s columns
    // TODO: (1) Optimize; and (2) try to avoid the if/else statement
    int u = Q_(rtemp - 1, rtemp - 1);
    Eigen::VectorXi a = A_(Eigen::all, rtemp - 1);
    A_.conservativeResize(Eigen::NoChange, rtemp - 1);
    if (rtemp == 1) { // Need to handle this special case separately because
                      // otherwise I have issues with Eigen...
        Q_.conservativeResize(0, 0);
    } else {
        // print();
        Eigen::VectorXi q = Q_.topRightCorner(rtemp - 1, 1);
        Q_.conservativeResize(rtemp - 1, rtemp - 1);
        Q_ += 2 * z * Eigen::MatrixXi(q.asDiagonal());
        ReduceQ();
    }
    b_ += z * a;
    ReduceVectorMod(b_, 2);
    phase_ = (phase_ + 2 * z * u) % 8;
    std::cout << "Almost done FixFinalBit";
    print();
    std::cout << "r = " << r_;
    pivots_.erase(rtemp - 1);
    --r_;
    std::cout << "Finished FixFinalBit";
    print();
}

void AffineState::ZeroColumnElim(int c) {
    // Recall: Matrix A_ has rank r_ and r_+1 columns at this stage
    std::cout << "Entering ZeroColumnElim...";
    print();
    ReindexSwapColumns(
        c, r_); // Step 1. We use r_ instead of r_+1 here because of 0-indexing
    std::cout << "After calling ReindexSwapColumns";
    print();
    Eigen::VectorXi q = Q_(Eigen::seq(0, r_ - 1), r_); // Step 2.
    ReduceVectorMod(q, 2); // TODO: Check that this reduction is appropriate.
    int u = Q_(r_, r_) % 4;
    A_.conservativeResize(Eigen::NoChange, r_);
    Q_.conservativeResize(r_, r_);
    pivots_.erase(r_); // TODO: IS THIS CORRECT???????
    if (u % 2 == 1) {
        Q_ += (u - 2) * q * q.transpose();
        ReduceQ();
        phase_ = (phase_ - u + 2) % 8;
        return;
    } else {
        std::cout << "Entering else part of ZeroColumnElim";
        print();
        std::cout << "q = " << q.transpose() << std::endl;
        // TODO: Check whether q == \vec{0}? For now assume not
        int ell;
        for (int i = 0; i < r_; i++) {
            if (q(i) != 0) {
                ell = i;
                break;
            }
        }
        std::cout << "Selected value of ell = " << ell << std::endl;
        for (int k = 0; k < r_; k++) {
            if (q(k) != 0 && k != ell) {
                ReindexSubtColumn(k, ell);
            }
        }
        std::cout << "Finished the ReindexSubtColumn part of ZeroColumnElim";
        print();
        ReindexSwapColumns(r_ - 1, ell); // r_-1 because of 0-indexing
        std::cout << "Finished the ReindexSwapColumns part of ZeroColumnElim";
        print();
        FixFinalBit(u / 2);
        // for (auto const& pair : pivots_) { // If column c was a pivot,
        //	if (pair.first == c) { pivots_.erase(pair.first); break; }
        // }
        return;
    }
}

// GATES:
void AffineState::H(int j) {
    ReduceQ(); // TODO: Removing this probably is fine
    std::cout << "State entering Hadamard subroutine:";
    print();
    // Step 1: Find c.
    int c = -1;
    for (int i = 0; i < r_; i++) {
        if (pivots_[i] == j) {
            c = i;
        }
    }
    std::cout << "Step 1: found c = " << c << std::endl;
    if (c > -1) { // Step 2. (Using -1 instead of 0 because of 0-indexing.)
        std::cout << "Step 2. Before ReselectPrincipalRow: ";
        print();
        ReselectPrincipalRow(j, c);
        std::cout << "Step 2. After ReselectPrincipalRow: ";
        print();
        if (j != pivots_[c]) {
            c = -1;
        }
        std::cout << "After Step 2:  c = " << c << std::endl;
    }

    // Step 3:
    std::cout << "Step 3. A = " << std::endl << A_ << std::endl;
    Eigen::VectorXi atilde;
    atilde.setZero(r_ + 1); // r_+1 entries
    atilde(Eigen::seq(0, r_ - 1)) =
        A_(j, Eigen::all)
            .transpose(); // Set the first r_ of them to row j of A_.
    std::cout << "Step 3. atilde = " << std::endl << atilde << std::endl;

    // Step 4:
    std::cout << "Before Step 4";
    print();
    A_.row(j).setZero();
    A_.conservativeResize(
        Eigen::NoChange,
        r_ + 1); // Add column, and (next two lines) set new column to e_j
    A_.col(r_).setZero();
    A_(j, r_) = 1;
    pivots_[r_] = j;
    std::cout << "After Step 4";
    print();

    // Step 5:
    std::cout << "Before Step 5";
    print();
    Q_.conservativeResize(r_ + 1, r_ + 1);
    Q_.row(r_) = atilde.transpose();
    Q_.col(r_) = atilde.transpose();
    Q_(r_, r_) = 2 * b_(j); // Already reduced mod 4
    std::cout << "After Step 5";
    print();

    // Step 6:
    b_(j) = 0;
    if (c > -1) {
        ZeroColumnElim(c);
    } else {
        ++r_;
    }
}

void AffineState::CZ(int j, int k) {
    if (A_.cols() > 0) { // If psi is a computational basis state then A_ has
                         // zero columns, so the next three lines do nothing
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
    Eigen::Vector<int, Eigen::Dynamic> a_j = A_(j, Eigen::all).transpose();
    Q_ += (1 - 2 * b_(j)) * a_j * a_j.transpose();
    ReduceQ();
    phase_ = (phase_ + 2 * b_(j)) % 8;
}

void AffineState::X(int j) { b_(j) = (b_(j) + 1) % 2; }

void AffineState::Z(int j) {
    phase_ = (phase_ + 4 * (b_(j))) % 8;
    Q_ = Q_ + 2 * Eigen::MatrixXi(A_(j, Eigen::all).asDiagonal());
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
            }
        }
        ReindexSwapColumns(k, r_ - 1);
        MakePrincipal(r_ - 1, j);
        FixFinalBit(beta);
        return beta;
    }
}

void ReduceMatrixMod(Eigen::MatrixXi& M,
                     int modulus) { // TODO: This is a really silly way of doing
                                    // it. Should find something better
    M = M.array() - (modulus * (M.array() / modulus));
    M = M.array() + modulus;
    M = M.array() - (modulus * (M.array() / modulus));
}

void ReduceVectorMod(Eigen::VectorXi& v,
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

std::ostream& operator<<(std::ostream& out, AffineState const& psi) {
    out << "n = " << psi.n_ << '\n';
    out << "r = " << psi.r_ << '\n';
    out << "phase = " << psi.phase_ << '\n';
    out << "Q = " << std::endl << psi.Q_ << '\n';
    out << "A = " << std::endl << psi.A_ << '\n';
    out << "b^T = " << psi.b_.transpose() << '\n';
    std::cout << "pivots = ";
    for (const auto& p : psi.pivots_) {
        out << "(" << p.first << ", " << p.second << "), ";
    }
    return out;
}

void AffineState::print() const { std::cout << '\n' << *this << std::endl; }

// int RandomBit() { // TODO: Optimize. It's inefficient to build a new
// generator each time a random bit is needed. 	std::random_device rd;
//	std::mt19937 mt(rd());
//	std::uniform_int_distribution<> distr(0, 1);
//	return distr(mt);
// }
