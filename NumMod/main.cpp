#include <iostream>
#include <fstream>

#include "Matrix.hpp"
#include "MatrixFunctions.hpp"

int main() {
    std::ifstream infile("matrix1.txt");
    sbl::Matrix<long double> m;
    infile >> m;

    //std::cout << sbl::det(m) << std::endl;

    //std::cout << m << std::endl << std::endl;

    /*sbl::GaussianElimination(m);
    std::cout << m << std::endl << std::endl;*/

    /*sbl::Matrix<long double>::LDU lu = sbl::LDUDecomposition(m);
    std::cout << lu.L << std::endl << std::endl;
    std::cout << lu.D << std::endl << std::endl;
    std::cout << lu.U << std::endl << std::endl;
    std::cout << lu.L * lu.D * lu.U << std::endl << std::endl;
    std::cout << sbl::det(m) << std::endl << std::endl;*/

    /*sbl::Inverse(m);
    std::cout << m << std::endl << std::endl;*/

    std::cout << m << std::endl << std::endl;
    //std::cout << HouseholderMatrix(m) << std::endl << std::endl;
    std::cout << norm1(m) << std::endl;
    std::cout << normFrob(m) << std::endl;
    std::cout << normInf(m) << std::endl;
    std::cout << (~m) * m << std::endl << std::endl;


    std::cin.ignore('\n', 1000);

    return 0;
}