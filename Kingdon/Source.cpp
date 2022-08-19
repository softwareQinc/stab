#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <stdlib.h>
#include "AffineState.h"

using namespace Eigen;

void aReduceGramRowCol(int c, MatrixXi & Q) {
	Q(c, c) = (4 + Q(c, c) % 4) % 4;
	for (int i = 0; i < Q.cols(); i++) {
		if (i != c) {
			Q(c, i) = (2 + Q(c, i) % 2) % 2;
			Q(i, c) = (2 + Q(i, c) % 2) % 2;
		}
	}
	return;
}

int main()
{
	/*MatrixXi Q;
	Q.setRandom(4, 4);
	std::cout << Q << std::endl;
	aReduceGramRowCol(3, Q);
	std::cout << Q << std::endl;*/


	AffineState psi(3);
	psi.X(0);
	psi.X(2);
	psi.H(0);
	psi.H(1);
	psi.CZ(0, 2);
	psi.H(0);
	psi.Z(2);
	psi.CZ(0, 1);
	/*psi.H(1);
	psi.H(2);
	psi.CZ(0, 1);
	psi.CZ(0, 2);
	psi.H(0);
	psi.H(1);
	psi.H(2);*/
	std::cout << psi;

	//for (int i = 1; i < psi.n; i++) {
	//	psi.CZ(0, i);
	//}

	//for (int i = 0; i < psi.n; i++) {
	//	psi.H(i);
	//}
	//for (int i = 1; i < psi.n; i++) {
	//	psi.CZ(0, i);
	//}
	//std::cout << psi;
	return 0;
}