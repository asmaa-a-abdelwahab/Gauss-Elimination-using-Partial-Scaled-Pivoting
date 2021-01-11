#include <iostream>
#include <vector>
#include "gauss_elimination.h"

using namespace std;

int main()
{
    int iDim;
    std::vector<float> _scale;
    std::cout.precision(4);        //set precision
    std::cout.setf(ios::fixed);    //display using the fixed precision
    std::cout << "\nEnter the no. of equations\n";
    std::cin >> iDim;                //input the no. of equations

    float *p[iDim];
    float arrA[iDim][iDim+1];  //Array declaration to store the augmented-matrix elements
    std::cout << "\nEnter the elements of the augmented-matrix row-wise:\n";
    for (int i = 0; i < iDim; i++)
        for (int j = 0; j <= iDim; j++)
        {
            std::cout << "A[" << i << "][" << j << "]=";
            std::cin >> arrA[i][j];
        };
    for (int i = 0; i < iDim; ++i)
        p[i] = arrA[i];

    GaussElimination matrix = GaussElimination();
    std::cout << "\nThis is the matrix to be solved\n";
    matrix.printMatrix(p, iDim);
    matrix.selectMax(p, iDim, _scale);
    for (int step = 0; step < iDim-1; step++)
        matrix.solve(p, iDim,_scale, step);
    matrix.backwardSubstitution(p, iDim);
    return 0;
}
