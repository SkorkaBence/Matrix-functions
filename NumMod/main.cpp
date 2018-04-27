#include <iostream>
#include <fstream>

#include "Matrix.hpp"
#include "MatrixFunctions.hpp"
#include "Polynomial.hpp"

int main() {
    std::ifstream infile("matrix1.txt");
    sbl::Matrix<long double> m;
    infile >> m;

    sbl::polynomial poly;
    poly.setValue(4, 1);
    poly.setValue(3, 1);
    poly.setValue(2, 5);
    poly.setValue(1, 3);
    poly.setValue(0, -10);

    std::cout << poly << std::endl << std::endl;

    sbl::Matrix<int> m2 = CompanionMatrix(poly);

    std::cout << m2 << std::endl << std::endl;

    std::cout << trace(m2) << std::endl << std::endl;

    std::cin.ignore('\n', 1000);

    return 0;
}