#include "AffineState.h"
#include <Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <vector>
#include <map>

using namespace Eigen;

using std::cout; // Just to make printing things during the debugging stage easier. Can remove later
using std::endl; // Just to make printing things during the debugging stage easier. Can remove later

// Initialization:
AffineState::AffineState(int m) {
	n = m;
	phase = 0;
	Q.setZero(0, 0);
	A.setZero(m, 0);
	b.setZero(m);
	r = 0;
	return;
}

// Subroutines:
void AffineState::ReduceGramRowCol(int c) { // TODO: Optimize
	Q(c, c) = (4 + Q(c, c) % 4) % 4;
	for (int i = 0; i < Q.cols(); i++) {
		if (i != c) {
			Q(c, i) = (2 + Q(c, i) % 2) % 2;
			Q(i, c) = (2 + Q(i, c) % 2) % 2;
		}
	}
	return;
}

void AffineState::ReindexSubtColumn(int k, int c) {
	if (c == k) { return; }
	else {
		for (int j = 0; j < n; j++) {
			if (A(j, c) != 0) { A(j, k) = (A(j, k) + 1) % 2; }
		}
		for (int h = 0; h < Q.cols(); h++) {
			if (Q(h, c) != 0) { Q(h, k) -= Q(h, c); }
			if (Q(c, h) != 0) { Q(k, h) -= Q(c, h); }
		}
		ReduceGramRowCol(k);
		//ReduceQ();
		return;
	}
}

void AffineState::ReindexSwapColumns(int k, int c) {
	if (c == k) { return; }
	else {
		int swaptmp;
		for (int j = 0; j < n; j++) {
			if (A(j, c) != 0 || A(j, k) != 0) {//TODO: Optimize swap
				swaptmp = A(j, c);
				A(j, c) = A(j, k);
				A(j, k) = swaptmp;
			}
		}
		for (int j = 0; j < Q.cols(); j++) {
			if (Q(j, c) != 0 || Q(j, k) != 0) {//TODO: Optimize swap
				swaptmp = Q(j, c);
				Q(j, c) = Q(j, k);
				Q(j, k) = swaptmp;
			}
		}
		for (int j=0; j < Q.cols(); j++){
			if (Q(c, j) != 0 || Q(k, j) != 0) {//TODO: Optimize swap
				swaptmp = Q(c, j);
				Q(c, j) = Q(k, j);
				Q(k, j) = swaptmp;
			}
		}
		ReduceQ();
		std::swap(pivots.at(k), pivots.at(c));
		return;
	}
}

void AffineState::MakePrincipal(int c, int j) {
	if (A(j, c) == 0) { std::cout << "Error......." << std::endl; return; }
	else {
		for (int k = 0; k < r; k++) {
			if (k != c && A(j, k) != 0) {
				ReindexSubtColumn(k, c);
			}
		}
		pivots[c] = j;
		return;
	}
}

void AffineState::ReselectPrincipalRow(int j, int c) {
	int j_star = -1; // Equivalent to j_star = 0 in the paper (but we need to use -1 because of 0-indexing)
	for (int jj = 0; jj < r; jj++) { // TODO: Rewrite to find optimal j_star (see paper)
		if (A(jj, c) != 0) { j_star = jj; break; }
	}
	if (j_star != -1) {
		MakePrincipal(c, j_star);
		return;
	}
	else { // TODO: Don't think is an error is necessary here, but should doublecheck eventually.
		return;
	}
}

void AffineState::FixFinalBit(int z) {
	cout << "Entering FixFinalBit"; print();
	int rtemp = A.cols(); // Define this because sometimes this function gets called from the ZeroColumnElim function, in which r might temporarily not be out of sync with A's columns
	// TODO: (1) Optimize; and (2) try to avoid the if/else statement
	int u = Q(rtemp - 1, rtemp - 1);
	VectorXi a = A(all, rtemp - 1);
	A.conservativeResize(NoChange, rtemp - 1);
	if (rtemp == 1) {// Need to handle this special case separately because otherwise I have issues with Eigen...
		Q.conservativeResize(0, 0);
	}
	else {
	//print();
		VectorXi q = Q.topRightCorner(rtemp - 1, 1);
		Q.conservativeResize(rtemp - 1, rtemp - 1);
		Q += 2 * z * MatrixXi(q.asDiagonal());
		ReduceQ();
	}
	b += z * a;
	ReduceVectorMod(b, 2);
	phase = (phase + 2 * z * u) % 8;
	cout << "Almost done FixFinalBit"; print();
	cout << "r = " << r;
	pivots.erase(rtemp - 1);
	--r;
	cout << "Finished FixFinalBit"; print();
	return;
}

void AffineState::ZeroColumnElim(int c) {
	// Recall: Matrix A has rank r and r+1 columns at this stage
	cout << "Entering ZeroColumnElim..."; print();
	ReindexSwapColumns(c, r); // Step 1. We use r instead of r+1 here because of 0-indexing
	cout << "After calling ReindexSwapColumns"; print();
	VectorXi q = Q(seq(0, r-1), r); // Step 2.
	ReduceVectorMod(q, 2); // TODO: Check that this reduction is appropriate.
	int u = Q(r, r) % 4;
	A.conservativeResize(NoChange, r);
	Q.conservativeResize(r, r);
	pivots.erase(r); // TODO: IS THIS CORRECT???????
	if (u % 2 == 1) {
		Q += (u - 2) * q * q.transpose();
		ReduceQ();
		phase = (phase - u + 2) % 8;
		return;
	}
	else {
		cout << "Entering else part of ZeroColumnElim"; print();
		cout << "q = " << q.transpose() << endl;
		// TODO: Check whether q == \vec{0}? For now assume not
		int ell;
		for (int i = 0; i < r; i++) {
			if (q(i) != 0) { ell = i; break; }
		}
		cout << "Selected value of ell = " << ell << endl;
		for (int k = 0; k < r; k++) {
			if (q(k) != 0 && k != ell) { ReindexSubtColumn(k, ell); }
		}
		cout << "Finished the ReindexSubtColumn part of ZeroColumnElim"; print();
		ReindexSwapColumns(r - 1, ell); // r-1 because of 0-indexing
		cout << "Finished the ReindexSwapColumns part of ZeroColumnElim"; print();
		FixFinalBit(u / 2);
		//for (auto const& pair : pivots) { // If column c was a pivot, 
		//	if (pair.first == c) { pivots.erase(pair.first); break; }
		//}
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
	for (int i = 0; i < r; i++) {
		if (pivots[i] == j) { c = i; }
	}
	std::cout << "Step 1: found c = " << c << std::endl;
	if (c > -1) { // Step 2. (Using -1 instead of 0 because of 0-indexing.)
		std::cout << "Step 2. Before ReselectPrincipalRow: "; print();
		ReselectPrincipalRow(j, c);
		std::cout << "Step 2. After ReselectPrincipalRow: "; print();
		if (j != pivots[c]) { c = -1; }
		std::cout << "After Step 2:  c = " << c << std::endl;
	}

	// Step 3:
	std::cout << "Step 3. A = " << std::endl << A << std::endl;
	VectorXi atilde;
	atilde.setZero(r + 1); // r+1 entries
	atilde(seq(0, r - 1)) = A(j, all).transpose(); // Set the first r of them to row j of A.
	std::cout << "Step 3. atilde = " << std::endl << atilde << std::endl;

	// Step 4:
	std::cout << "Before Step 4"; print();
	A.row(j).setZero();
	A.conservativeResize(NoChange, r + 1); // Add column, and (next two lines) set new column to e_j
	A.col(r).setZero();
	A(j, r) = 1;
	pivots[r] = j;
	std::cout << "After Step 4"; print();

	// Step 5:
	std::cout << "Before Step 5"; print();
	Q.conservativeResize(r + 1, r + 1);
	Q.row(r) = atilde.transpose();
	Q.col(r) = atilde.transpose();
	Q(r, r) = 2 * b(j); // Already reduced mod 4
	std::cout << "After Step 5"; print();

	// Step 6:
	b(j) = 0;
	if (c > -1) { ZeroColumnElim(c); }
	else { ++r; }
	return;
}

void AffineState::CZ(int j, int k) {
	if (A.cols() > 0) { // If psi is a computational basis state then A has zero columns, so the next three lines do nothing
		VectorXi a_j = A(j, all).transpose();
		VectorXi a_k = A(k, all).transpose();
		Q += a_j * a_k.transpose() + a_k * a_j.transpose() + 2 * MatrixXi((b(k) * a_j.asDiagonal() + b(j) * a_k.asDiagonal()));
		ReduceQ();
	}
	phase = (phase + 4 * b(j) * b(k)) % 8;
	return;
}

void AffineState::CX(int j, int k) { // TODO: Use non-naive method.
	H(k);
	CZ(j, k);
	H(k);
	return;
}

void AffineState::S(int j) {
	Vector<int, Dynamic> a_j = A(j, all).transpose();
	Q += (1 - 2 * b(j)) * a_j * a_j.transpose();
	ReduceQ();
	phase = (phase + 2 * b(j)) % 8;
	return;
}

void AffineState::X(int j) {
	b(j) = (b(j) + 1) % 2;
	return;
}

void AffineState::Z(int j) {
	phase = (phase + 4 * (b(j))) % 8;
	Q = Q + 2 * MatrixXi(A(j, all).asDiagonal());
	ReduceQ();
	return;
}

void AffineState::Y(int j) {
	phase = (phase + 2) % 8;
	Z(j);
	X(j);
	return;
}

int AffineState::MeasureZ(int j) {
	// Check whether deterministic:
	if (A(j, all).isZero()) {
		return b(j);
	}
	else {
		int beta = std::rand() % 2; // TODO: Optimize and use better source of randomness
		int k;
		for (int i = 0; i < r; i++) { // TODO: Optimize
			if (A(j, i) != 0) { k = i; }
		}
		ReindexSwapColumns(k, r - 1);
		MakePrincipal(r - 1, j);
		FixFinalBit(beta);
		return beta;
	}
}

void ReduceMatrixMod(MatrixXi& M, int modulus) { // TODO: This is a really silly way of doing it. Should find something better
	M = M.array() - (modulus * (M.array() / modulus));
	M = M.array() + modulus;
	M = M.array() - (modulus * (M.array() / modulus));
	return;
}

void ReduceVectorMod(VectorXi& v, int modulus) { // TODO: This is a really silly way of doing it. Should find something better
	v = v.array() - (modulus * (v.array() / modulus));
	v = v.array() + modulus;
	v = v.array() - (modulus * (v.array() / modulus));
	return;
}

void AffineState::ReduceQ() { // TODO: Optimize
	// Need to reduce off-diagonal elements modulo 2 and diagonal elements modulo 4
	/*std::cout << "Reducing Q..." << std::endl << std::endl;
	std::cout << "Q = " << std::endl << Q << std::endl;*/
	MatrixXi Qdiag = Q.diagonal().asDiagonal();
	Q = Q - Qdiag;
	ReduceMatrixMod(Q, 2);
	ReduceMatrixMod(Qdiag, 4);
	Q += Qdiag;
	/*std::cout << "Now, Q = " << std::endl << Q << std::endl;*/
	return;
}

std::ostream& operator<< (std::ostream& out, AffineState const& psi) {
	out << "n = " << psi.n << std::endl;
	out << "r = " << psi.r << std::endl;
	out << "phase = " << psi.phase << std::endl;
	out << "Q = " << std::endl << psi.Q << std::endl;
	out << "A = " << std::endl << psi.A << std::endl;
	out << "b^T = " << psi.b.transpose() << std::endl;
	std::cout << "pivots = ";
	for (const auto& p : psi.pivots) {
		std::cout << "(" << p.first << ", " << p.second << "), ";
	}
	return out;
}

void AffineState::print() {
	std::cout << std::endl;
	std::cout << "n = " << n << std::endl;
	std::cout << "r = " << r << std::endl;
	std::cout << "phase = " << phase << std::endl;
	std::cout << "Q = " << std::endl << Q << std::endl;
	std::cout << "A = " << std::endl << A << std::endl;
	std::cout << "b^T = " << b.transpose() << std::endl;
	std::cout << "pivots = ";
	for (const auto& p : pivots) {
		std::cout << "(" << p.first << ", " << p.second << "), ";
	}
	std::cout << std::endl;
	return;
}

//int RandomBit() { // TODO: Optimize. It's inefficient to build a new generator each time a random bit is needed.
//	std::random_device rd;
//	std::mt19937 mt(rd());
//	std::uniform_int_distribution<> distr(0, 1);
//	return distr(mt);
//}
