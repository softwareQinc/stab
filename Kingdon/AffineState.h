#pragma once
#include <Eigen/Dense>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <random>

using namespace Eigen;

class AffineState
{
public: 
	// State is represented by exp(i*pi*phase/8)/sqrt(2^r) \sum_{x \in \{0,1\}^r} i^{x^T Q x} \ket{Ax + b mod 2}, where A is n\times r and has rank r <= n.
	int n; // number of qubits
	int phase;  // global phase is exp(i*pi*phase/8)
	Matrix<int, Dynamic, Dynamic> Q; // quadratic function
	Matrix<int, Dynamic, Dynamic> A; // affine space generating matrix
	Vector<int, Dynamic> b; // affine space offset
	std::map<int, int> pivots; // AKA "principal index map." Keys are columns, values are the rows that contain pivots in those columns. Note that we are using zero-indexing, so the smallest key (assuming the map is nonempty) will always be 0, and the corresponding value is the index of the row that has a pivot in column 0.
	

	AffineState(int m);


	// Gates and measurements
	void CZ(int a, int b);
	void CX(int a, int b);
	void S(int j);
	void H(int a);
	void X(int j);
	void Y(int j);
	void Z(int j);
	int MeasureZ(int j); // Returns outcome and updates state


	// Option to print state
	friend std::ostream& operator<< (std::ostream& out, AffineState const& psi);

private:
	int r; // Technically unnecessary since this is just A.cols(), but it is handy to not have to declare it each time


	// Subroutines
	void FixFinalBit(int z);
	void ReduceGramRowCol(int c);
	void ReindexSubtColumn(int k, int c);
	void ReindexSwapColumns(int k, int c);
	void MakePrincipal(int c, int j);
	void ReselectPrincipalRow(int j, int c);
	void ZeroColumnElim(int c);


	void ReduceQ(); // Reduces mod 4 on the diagonal and mod 2 on the off-diagonal


	void print(); // Only used for printing state at intermediate stages of the calculation when diagnosing errors
};

// Small helper functions. Not sure if this is the best place to declare them?
void ReduceMatrixMod(MatrixXi& M, int modulus);
void ReduceVectorMod(VectorXi& v, int modulus);
