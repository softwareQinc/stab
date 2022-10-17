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
        Qmaster_.setZero(n + 1, n + 1);
        Qonesmaster_.setZero(n + 1);
        Qtwosmaster_.setZero(n + 1);
        Amaster_.setZero(n, n + 1); // We sometimes need an extra column for workspace
        Q_ = std::make_unique<block_t>(Qmaster_, 0, 0, 0, 0);
        Qones_ = std::make_unique<subvec_t>(Qonesmaster_, 0, 0);
        Qtwos_ = std::make_unique<subvec_t>(Qtwosmaster_, 0, 0);
        A_ = std::make_unique<block_t>(Amaster_, 0, 0, n, 0);
        b_.setZero(n);
    }

    // copy
    AffineState::AffineState(const AffineState &other) { // TODO: How to handle pointers?
        n_ = other.n_;
        phase_ = other.phase_;
        r_ = other.r_;
        Qmaster_ = other.Qmaster_;
        Qonesmaster_ = other.Qonesmaster_;
        Qtwosmaster_ = other.Qtwosmaster_;
        Amaster_ = other.Amaster_;
        b_ = other.b_;
        pivots_ = other.pivots_;

        //TODO check this
        Q_ = std::make_unique<block_t>(Qmaster_, 0, 0, r_, r_);
        Qones_ = std::make_unique<subvec_t>(Qonesmaster_, 0, r_);
        Qtwos_ = std::make_unique<subvec_t>(Qtwosmaster_, 0, r_);
        A_ = std::make_unique<block_t>(Amaster_, 0, 0, n_, r_);
    }

    int AffineState::n() const { return n_; }

    int AffineState::phase() const { return phase_; }

    mat_b_t AffineState::Q() const { return *Q_; }// TODO: Fix these

    mat_b_t AffineState::A() const { return *A_; }

    vec_b_t AffineState::b() const { return b_; }

    std::unordered_map<int, int> AffineState::pivots() const { return pivots_; }

    int AffineState::r() const { return r_; }

// Subroutines:
    void AffineState::DiagAdd(int index, int addend) { // TODO: Doublecheck
        addend %= 4; // TODO: Might be able to remove this
        bool ones = bool(addend % 2);
        bool twos = (addend > 1);

        (*Qtwos_)(index) ^= (twos ^ ((*Qones_)(index)&ones));
        (*Qones_)(index) ^= ones;
    }

    void AffineState::ReindexSubtColumn(int k, int c) {
        // Applies the update Column k of A <--- Column k + Column c
        assert(c != k);
        (*Qtwos_)(k) ^= (*Q_)(k, c);
        DiagAdd(k, (4 + 2 * int((*Qtwos_)(c)) + int((*Qones_)(c)) % 4)); // TODO: Doublecheck the sign
        
        A_->col(k) = A_->col(k).binaryExpr(A_->col(c), bool_xor);
        Q_->col(k) = Q_->col(k).binaryExpr(Q_->col(c), bool_xor);
        Q_->row(k) = Q_->row(k).binaryExpr(Q_->row(c), bool_xor);
        (*Q_)(c, k) ^= (*Qones_)(c);
        (*Q_)(k, c) ^= (*Qones_)(c);

        (*Q_)(k, k) = false;
    }

    void AffineState::ReindexSwapColumns(int k, int c) {
        if (c != k) {
            A_->col(k).swap(A_->col(c));
            Q_->col(k).swap(Q_->col(c));
            Q_->row(k).swap(Q_->row(c));
            std::swap((*Qones_)(k), (*Qones_)(c));
            std::swap((*Qtwos_)(k), (*Qtwos_)(c));
            std::swap(pivots_.at(k), pivots_.at(c));
        }
    }

    void AffineState::MakePrincipal(int c, int j) {
        assert((*A_)(j, c));
        for (int k = 0; k < r_; ++k) {
            if ((*A_)(j, k) && k != c) {
                ReindexSubtColumn(k, c);
            }
        }
        pivots_[c] = j;
    }

    bool AffineState::ReselectPrincipalRow(int j, int c) {
        int j_star = -1; // Equivalent to j_star = 0 in the paper (but we need to
        // use -1 because of 0-indexing)
        for (int jj = 0; jj < n_; ++jj) {
            if ((*A_)(jj, c) && jj != j) {
                j_star = jj;
                break;
            }
        }

        if (j_star != -1) {
            MakePrincipal(c, j_star);
            return true;
        } else {
            return false;
        }
    }

    void AffineState::FixFinalBit(int z) {
        assert(r_ != 0);

        // Step 1:
        vec_b_t a = A_->col(r_ - 1); // r_ - 1 because of 0-indexing
        vec_b_t q;
        if (r_ > 1) { // If r == 1 then the line below throws an error due to r_ - 2
            q = (*Q_)(Eigen::seq(0, r_ - 2), r_ - 1);
        }
        int u = 2 * int((*Qtwos_)(r_ - 1)) + int((*Qones_)(r_ - 1));

        // Step 2:
        A_->col(r_ - 1).setZero(); // Set to zero before removing since this col will persist in Amaster_
        Q_->col(r_ - 1).setZero(); // Ditto
        Q_->row(r_ - 1).setZero(); // Ditto
        (*Qtwos_)(r_ - 1) = false; // Ditto
        (*Qones_)(r_ - 1) = false; // Ditto
        A_ = std::make_unique<block_t>(Amaster_, 0, 0, n_, r_ - 1);
        Q_ = std::make_unique<block_t>(Qmaster_, 0, 0, r_ - 1, r_ - 1);
        Qones_ = std::make_unique<subvec_t>(Qonesmaster_, 0, r_ - 1);
        Qtwos_ = std::make_unique<subvec_t>(Qtwosmaster_, 0, r_ - 1);

        // Step 3:
        if (z == 1) { // ReduceMod is relatively expensive, so it makes sense to check whether z == 1
            for (int i = 0; i < r_ - 1; ++i) {
                if (q(i)) {
                    (*Qtwos_)(i) ^= 1;
                }
            }
            b_ = b_.binaryExpr(a, bool_xor);
            phase_ = (phase_ + 2 * u) % 8;
        }

        // Step 4:
        pivots_.erase(r_ - 1);
        --r_;
    }

    void AffineState::ZeroColumnElim(int c) {
        // Recall: Matrix A_ has rank r_ and r_+1 columns at this stage

        // Step 1:
        ReindexSwapColumns(c, r_); // Use r_ instead of r_+1 here because of 0-indexing

        // Step 2:
        vec_b_t q = (*Q_)(Eigen::seq(0, r_ - 1), r_); // Step 2.
        int u = 2 * int((*Qtwos_)(r_)) + int((*Qones_)(r_));

        // Step 3:
        // No need to reset A_->col(r_) to zero since it's already zero
        A_ = std::make_unique<block_t>(Amaster_, 0, 0, n_, r_);
        Q_->row(r_).setZero();
        Q_->col(r_).setZero();
        (*Qones_)(r_) = false;
        (*Qtwos_)(r_) = false;
        Q_ = std::make_unique<block_t>(Qmaster_, 0, 0, r_, r_);
        Qones_ = std::make_unique<subvec_t>(Qonesmaster_, 0, r_);
        Qtwos_ = std::make_unique<subvec_t>(Qtwosmaster_, 0, r_);

        // Step 3.5:
        pivots_.erase(r_); // <-- Both parts of the if statement in Step 4 require us to do this

        // Step 4
        if (u % 2 == 1) {
            (*Q_) = (*Q_).binaryExpr(q * q.transpose(), bool_xor);
            Q_->diagonal().setZero();
            for (int i = 0; i < r_; ++i) {
                if (q(i)) {
                    DiagAdd(i, u + 2);
                }
            }
            phase_ = (phase_ - u + 2) % 8;
            return;
        } else { // u is even
            assert(!q.isZero());
            int ell;
            for (int i = 0; i < r_; ++i) {
                if (q(i)) {
                    ell = i;
                    break;
                }
            }
            for (int k = 0; k < r_; ++k) {
                if (q(k) && k != ell) {
                    ReindexSubtColumn(k, ell);
                }
            }
            ReindexSwapColumns(r_ - 1, ell); // r_-1 because of 0-indexing. This step also produces A^(3) and Q^(3)
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
        vec_b_t atilde = vec_b_t::Zero(r_ + 1);
        atilde(Eigen::seq(0, r_ - 1)) = A_->row(j).transpose();
        vec_b_t atildetrans = atilde.transpose();

        // Step 4:
        A_->row(j).setZero();
        A_ = std::make_unique<block_t>(Amaster_, 0, 0, n_, r_ + 1);
        // Note that the new column of A_ is already zero
        (*A_)(j, r_) = true;
        pivots_[r_] = j;

        // Step 5:
        Q_ = std::make_unique<block_t>(Qmaster_, 0, 0, r_ + 1, r_ + 1);
        Q_->row(r_) = atildetrans;
        Q_->col(r_) = atildetrans;
        Qones_ = std::make_unique<subvec_t>(Qonesmaster_, 0, r_ + 1);
        Qtwos_ = std::make_unique<subvec_t>(Qtwosmaster_, 0, r_ + 1);
        DiagAdd(r_, 2*int(b_(j)));

        // Step 6:
        b_(j) = false;

        // Step 7:
        if (c > -1) { // Case 2 from Alex's notes
            ZeroColumnElim(c);
        } else { // Case 1 from Alex's notes
            ++r_;
        }
    }

    void AffineState::CZ(int j, int k) {
        if (r_ > 0) { // If psi is a computational basis state then A_ has
            // zero columns, so this block would do nothing anyway

            int pcj = piv_col(j);
            int pck = piv_col(k);

            if (pcj == -1 && pck != -1) {
                std::swap(pcj, pck);
                std::swap(j, k);
            }

            if (pcj != -1 && pck != -1) { // Nice case since it takes time O(1)
                (*Q_)(pcj, pck) ^= 1;
                (*Q_)(pck, pcj) ^= 1;
                (*Qtwos_)(pcj) ^= b_(k);
                (*Qtwos_)(pck) ^= b_(j);

            } else if (pcj != -1 && pck == -1) { // O(r_) time
                Q_->col(pcj) =
                    Q_->col(pcj).binaryExpr(A_->row(k).transpose(), bool_xor);
                Q_->row(pcj) = Q_->row(pcj).binaryExpr(A_->row(k), bool_xor);
                (*Qtwos_)(pcj) ^= b_(k);
                for (int i = 0; i < r_; ++i) {
                    (*Qtwos_)(i) ^= b_(j) & (*A_)(k, i);
                }

            } else { // Slow case since it takes time O(r_^2)
                (*Q_) = (*Q_).binaryExpr(A_->row(j).transpose() * A_->row(k),
                                         bool_xor);
                (*Q_) = (*Q_).binaryExpr(A_->row(k).transpose() * A_->row(j),
                                         bool_xor);
                Q_->diagonal().setZero();
                for (int i = 0; i < r_; ++i) {
                    (*Qtwos_)(i) ^= (*A_)(j, i) & (*A_)(k, i);
                    (*Qtwos_)(i) ^= b_(k) & (*A_)(j, i);
                    (*Qtwos_)(i) ^= b_(j) & (*A_)(k, i);
                }
            }
        }
        phase_ = (phase_ + 4 * (b_(j) & b_(k))) % 8;
    }

    void AffineState::CX(int h, int j) {
        // Step 1:
        int c = piv_col(j);

        // Step 2:
        A_->row(j) = A_->row(j).binaryExpr(A_->row(h), bool_xor);

        // Step 3:
        b_(j) ^= b_(h);

        // Step 4:
        if (c != -1) {
            ReselectPrincipalRow(-1, c);
        }
    }

    void AffineState::SWAP(int j, int k) {
        // TODO: Probably better way to update pivots_
        A_->row(j).swap(A_->row(k));
        b_.row(j).swap(b_.row(k));
        int cskip = -1;
        for (int c = 0; c < r_; ++c) {
            if (pivots_[c] == j) {
                pivots_[c] = k;
                cskip = c;
                break;
            }
        }
        for (int c = 0; c < r_; ++c) {
            if (pivots_[c] == k && c != cskip) {
                pivots_[c] = j;
                break;
            }
        }
    }

    void AffineState::S_or_SDG(int j, bool dg) {
        // dg = true means we apply S^\dagger, dg = false means we apply S
        int sign = 1 - 2 * int(dg);  // = +1 for S and -1 for S^\dagger

        if (r_ > 0) { // The lines below are irrelevant if r_ == 0
            int pivcol = piv_col(j);
            if (pivcol != -1) { // Special case that can be done quickly
                DiagAdd(pivcol, (4 + sign * (1 - 2 * int(b_(j)))) % 4);
            } else { // In general takes O(r_) time

                (*Q_) = (*Q_).binaryExpr(A_->row(j) * A_->row(j).transpose(),
                                         bool_xor);
                for (int i = 0; i < r_; ++i) {
                    if ((*Q_)(i, i)) {
                        DiagAdd(i, 4 + sign * (1 - 2 * b_(j)));
                        (*Q_)(i, i) = false;
                    }
                }
            }
        }
        phase_ = (phase_ + sign * 2 * int(b_(j)) + 8) % 8;
    }

    void AffineState::S(int j) { S_or_SDG(j, false); }

    void AffineState::SDG(int j) { S_or_SDG(j, true); }

    void AffineState::X(int j) { b_(j) ^= 1; }

    void AffineState::Z(int j) {
        phase_ = (phase_ + 4 * int(b_(j))) % 8;
        *Qtwos_ = (*Qtwos_).binaryExpr(A_->row(j).transpose(), bool_xor);
    }

    void AffineState::Y(int j) {
        phase_ = (phase_ + 2) % 8;
        Z(j);
        X(j);
    }

    bool AffineState::MeasureZ(int j, bool postselect,
                              bool postselected_outcome) {
        // Default parameters are postselect=false, postselected_outcome=bool
        if (A_->row(j).isZero()) { // Deterministic case
            if (postselect && b_(j) != postselected_outcome) {
                throw std::logic_error("Postselection impossible");
            } else {
                return postselected_outcome;
            }
        } else { // Random case
            bool beta;
            if (postselect) {
                beta = postselected_outcome;
            } else {
                beta = bool(random_bit());
            }
            int k;
            for (int i = 0; i < r_; i++) { // TODO: Optimize
                if ((*A_)(j, i)) {
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

    std::vector<int> AffineState::MeasureAll() const {
        vec_b_t x(r_);
        // x.setZero(r_); // no need, since we sized it first
        std::for_each(x.data(), x.data() + x.size(), [](auto &elem) {
            elem = random_bit();
        });
        vec_b_t tmp = b_;
        // Matrix multiplication with vec_b_t doesn't use XOR, so do it "manually"
        for (int row = 0; row < n_; ++row) {
            for (int col = 0; col < r_; ++col) {
                tmp[row] ^= (*A_)(row, col) & x(col);
            }
        }

        return {tmp.data(), tmp.data() + tmp.size()}; // TODO: Doublecheck this
    }

    std::map<std::vector<int>, int>
    AffineState::Sample(int nreps) const { // TODO: Naive approach can be improved
        std::map<std::vector<int>, int> results;
        for (int repnumber = 0; repnumber < nreps; ++repnumber) {
            ++results[MeasureAll()];
        }
        return results;
    }

    void AffineState::Reset(int j) {
        bool res = MeasureZ(j);
        if (res) {
            X(j);
        }
    }

    Eigen::VectorXcd AffineState::to_vec() const {
        if (n_ > 16) {
            throw std::logic_error("Maximum number of qubits for statevector "
                                   "representation is 15");
        }


        Eigen::VectorXcd vec(int(pow(2, n_))); // This will be the statevector
        for (int k = 0; k < pow(2, r_); ++k) {
            // For each x \in \{0,1\}^ncols, we now add the corresponding term to vec
            // First, let x = vector whose entries are the binary digits of k:
            vec_b_t x;
            x.setZero(r_);
            std::bitset<16> bs(k);             // Get binary string
            for (int j = 0; j < r_; ++j) { // Cast binary string to x
                x(j) = bs[j];
            }

            // Figure out which basis state x results in:
            vec_b_t ket = b_;
            for (int row = 0; row < n_; ++row) {
                for (int col = 0; col < r_; ++col) {
                    ket(row) ^= (*A_)(row, col) & x(col);
                }
            }
            // "ket" is the binary representation of some number basis_state_number. We need to add the correct amplitude into vec[basis_state_number]. Note that AffineState puts the zero-th qubit on the left, but most people put it on the right, so we also make this conversion.

            int basis_state_number = 0;
            for (int j = 0; j < n_; ++j) {
                basis_state_number += int(int(ket(j)) * pow(2, n_ - 1 - j)); // Conversion
            }

            std::complex<double> phase = 0;
            for (int row = 0; row < r_; ++row) {
                for (int col = 0; col < r_; ++col) {
                    phase += 2 * int((*Q_)(row, col) & x(row) & x(col));
                }
                phase += 2 * int((*Qtwos_)(row)) + int((*Qones_)(row));
            }

            std::complex<double> pi = std::atan(1.0) * 4.0;
            std::complex<double> i(0, 1); // TODO: Better way of getting i and pi?
            vec[basis_state_number] += std::exp(i * pi * phase / 4.0);
        }
        return vec / pow(2, r_ / 2.0);
    }

    std::ostream &operator<<(std::ostream &out, AffineState const &psi) {
        out << "QUADRATIC FORM REPRESENTATION: \n";
        out << "phase = exp(" << psi.phase_ << "*i*pi/4)\n";
        out << "r = " << psi.r_ << "\n";

        Eigen::MatrixXi Q, A;
        Eigen::VectorXi b(psi.n_);
        Q.setZero(psi.r_, psi.r_);
        A.setZero(psi.n_, psi.r_);
        for (int row = 0; row < psi.r_; ++row) {
            for (int col = 0; col < psi.r_; ++col) {
                Q(row, col) = int((*psi.Q_)(row, col));
            }
            Q(row, row) =
                2 * int(psi.Qtwosmaster_(row)) + int(psi.Qonesmaster_(row));
        }
        for (int row = 0; row < psi.n_; ++row) {
            for (int col = 0; col < psi.r_; ++col) {
                A(row, col) = psi.Amaster_(row, col);
            }
            b(row) = int(psi.b_(row));
        }


        if (psi.r_ > 0) { // If psi is a computational basis state, no need to print this stuff
            out << "Q = " << std::endl << Q << '\n';
            out << "A = " << std::endl << A << '\n';
        }
        out << "b^T = " << b.transpose() << '\n';
        out << "pivots = ";
        for (auto p: psi.pivots_) {
            out << "(" << p.first << ", " << p.second << "), ";
        }

        //std::cout << "\n\n BASIS STATE DECOMPOSITION:\n";

        //// Very naive way of printing the state, but this will mainly be used
        //// for debugging. We can always optimize it later if it's important.
        //if (psi.r_ > 12) {
        //    std::cout << "Too many amplitudes to print\n";
        //} else {
        //    vec_u_t x;
        //    x.setZero(psi.r_);
        //    for (int i = 0; i < pow(2, psi.r_); ++i) {
        //        std::bitset<16> bs(i);         // Get binary string
        //        for (int j = 0; j < psi.r_; ++j) { // Cast binary string to x
        //            x(j) = int(bs[j]);
        //        }

        //        vec_u_t ket = (*psi.A_) * x + psi.b_;
        //        ket = ReduceMod(ket, 2);
        //        std::string rel_phase;
        //        int ampl = (2 * (x.transpose() * (*psi.Q_) * x)[0] + psi.phase_) % 8;
        //        out << "exp(" << ampl << "i*pi/4) * "
        //            << "|" << ket.transpose() << ">\n";
        //    }
        //}
        return out;
    }

    void AffineState::print() const { std::cout << '\n' << *this << std::endl; }

} // namespace stab
