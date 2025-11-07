#include <iostream>
#include <iomanip>
#include <cmath>

int const N = 5;

double A[N][N] = {
    {5, 4, 3, 2, 1},
    {10, 8, 7, 6, 5},
    {-1, 2, -3, 4, -5},
    {6, 5, -4, 3, -2},
    {1, 2, 3, 4, 5}
};

double b[N] = {37, 99, -9, 12, 53};
double x[N];
int rowsIndexes[N] = {0, 1, 2, 3, 4};

//********************************************************
//* Prints current A matrix with applied row permutations.
//********************************************************
void print_matrix_A()
{
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Current A matrix:\n";
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            std::cout << std::setw(10) << A[rowsIndexes[i]][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n=======================================\n\n";
}

//********************************************************
//* Performs LU decomposition of A with partial pivoting.
//********************************************************
void LU_decompose()
{
    std::cout << "LU decomposition:\n";

    for(int j = 0; j < N; j++)
    {
        // Pivoting: find max element in column if diagonal element is zero
        if(A[rowsIndexes[j]][j] == 0.0)
        {
            int maxRow = j;
            double max = fabs(A[j][j]);

            for(int i = j + 1; i < N; i++)
            {
                double a = fabs(A[i][j]);
                if(a > max)
                {
                    max = a;
                    maxRow = i;
                }
            }

            // Swap row indexes to perform partial pivoting
            int tmp = rowsIndexes[j];
            rowsIndexes[j] = rowsIndexes[maxRow];
            rowsIndexes[maxRow] = tmp;
        }

        // Gaussian elimination
        for(int w = j + 1; w < N; w++)
        {
            // Multiplier
            double wsp = A[rowsIndexes[w]][j] / A[rowsIndexes[j]][j];

            // Store L coefficient in A
            A[rowsIndexes[w]][j] = wsp;

            // Update remaining elements of the row
            for(int i = j + 1; i < N; i++)
            {
                A[rowsIndexes[w]][i] -= A[rowsIndexes[j]][i] * wsp;
            }
        }

        std::cout << "Matrix A after column elimination: " << j << ":\n";
        print_matrix_A();
    }
}

//**********************************************************
//* Solves Ly = b and Ux = y using the LU-decomposed matrix.
//**********************************************************
void LU_solve()
{
    double y[5] = {0, 0, 0, 0, 0};

    // Forward substitution: solve Ly = b
    for(int i = 0; i < N; i++)
    {
        y[i] = b[rowsIndexes[i]];
        for(int j = 0; j < i; j++)
        {
            y[i] -= A[rowsIndexes[i]][j] * y[j];
        }
    }

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\ny after solving Ly = b\n";
    for(int i = 0; i < N; ++i)
    {
        std::cout << "y[" << rowsIndexes[i] << "] = " << y[rowsIndexes[i]] << std::endl;
    }

    // Back substitution: solve Ux = y
    for(int i = 4; i >= 0; i--)
    {
        x[i] = y[i];
        for(int j = i + 1; j < N; j++)
        {
            x[i] -= A[rowsIndexes[i]][j] * x[j];
        }
        x[i] /= A[rowsIndexes[i]][i];
    }
}

int main()
{
    print_matrix_A();

    LU_decompose();
    LU_solve();

    std::cout << "\nAx = b solution:\n";
    for (int i = 0; i < N; ++i)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    return 0;
}