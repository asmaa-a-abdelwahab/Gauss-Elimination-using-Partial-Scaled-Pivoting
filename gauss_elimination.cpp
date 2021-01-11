//Gauss Elimination
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include "gauss_elimination.h"

using namespace std;

void GaussElimination::printMatrix(float **a, int n)
{
    int i,j;
    for (i=0;i<n;i++)            //print the new matrix
    {
        for (j=0;j<=n;j++)
            std::cout << a[i][j] << std::setw(16);
        std::cout << "\n";
    };
    std::cout << "\n";
}

static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}

void GaussElimination::selectMax(float **a, int n, vector<float> &scale)
{
    for (int i = 0; i < n; i++)
    {
        float max = *std::max_element(a[i], a[i] + n, abs_compare);
        scale.push_back(max);
    }
}

void GaussElimination::scaledPartialPivoting(float **a, int n, vector<float> &scale, int step)
{
    std::vector<float> scaledElements;
    for (int i = step; i < n; i++)
        scaledElements.push_back(std::abs(a[i][step]/scale[i]));
        //std::cout << a[i][step] << "\t/ scale\t" << scale[i] << "\t= scaledElement  (" << scaledElements[i-step] << ")\n" ;
    int maxScaledElementidx = std::max_element(scaledElements.begin(),scaledElements.end()) - scaledElements.begin();
    //std::cout << "\n" << maxScaledElementidx+step << "\n\n";
    if (maxScaledElementidx+step > step)
    {
        std::swap(a[step], a[maxScaledElementidx+step]);
        std::swap(scale[step], scale[maxScaledElementidx+step]);
    }
    //std::printf("This is the matrix after selecting num %d pivotEQ equation\n", step+1);
    //printMatrix(a, n);
}

void GaussElimination::forwardElimination(float **a, int n, int step)
{
    for (int i = step+1; i < n; i++)
    {
        float x = a[i][step];
        for (int j = step; j <= n; j++)
            a[i][j] = a[i][j] - ((x/a[step][step]) * a[step][j]);
    }

    std::printf("\nThe matrix after step %d of elimination\n", step+1);
    printMatrix(a, n);
}

void GaussElimination::solve(float **a, int n, vector<float> &scale, int step)
{
    scaledPartialPivoting(a, n, scale, step);
    forwardElimination(a, n, step);
}

void GaussElimination::backwardSubstitution(float **a, int n)
{
    float x[n], R[n];
    for (int i = n-1; i >= 0; i--)
    {
        x[i] = a[i][n];
        for (int j = i+1; j < n; j++)
            if (j != i)
                x[i] = x[i] - a[i][j] * x[j];
        x[i] = x[i] / a[i][i];
        std::printf("X%d =\t %.4f\n", i+1, x[i]);
   }
   std::cout << "\n";
   //R = AX-B
   for (int i = 0; i < n; i++)
       for (int j = 0; j < n; j++)
           R[i]+=( a[i][j]*x[j]);

   for (int i = 0; i < n; i++)
   {
       R[i] = R[i] - a[i][n];
       std::printf("R%d =\t %.4f\n", i+1, R[i]);
   }
   //check solution validity
   if (*max_element(R , R + n) <= 0.04)
       std::cout << "\nvalid_solution = true\n\n";
   else
       std::cout << "\nvalid_solution = false\n\n";
}
