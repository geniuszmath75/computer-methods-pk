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
//*
//* Funkcja wyświetlająca bieżący stan macierzy A z uwzględnieniem permutacji
//* wierszy.
//* Wykorzystuje tablicę rowsIndexes do wyświetlenia wierszy w odpowiedniej
//* kolejności
//* po wykonanych operacjach częściowego wyboru elementu podstawowego (partial pivoting).
//********************************************************

void print_matrix_A()
{
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Current A matrix:\n";
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            // Wyświetlenie elementu A[rowsIndexes[i]][j] - dostęp do wiersza
            // odbywa się przez tablicę indeksów, co umożliwia logiczną permutację
            // bez fizycznego przesuwania danych w pamięci
            std::cout << std::setw(10) << A[rowsIndexes[i]][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n=======================================\n\n";
}

//********************************************************
//* Performs LU decomposition of A with partial pivoting.
//*
//* Procedura realizująca dekompozycję LU macierzy A metodą eliminacji Gaussa
//* z częściowym wyborem elementu podstawowego (partial pivoting).
//* Wynik: macierz A zostaje zastąpiona przez macierze L i U zapisane w jednej strukturze:
//* - elementy na przekątnej i powyżej to macierz górnotrójkątna U
//* - elementy poniżej przekątnej to mnożniki (multipliers) macierzy dolnotrójkątnej L
//* - przekątna macierzy L jest jedynkowa (nie jest przechowywana jawnie)
//********************************************************
void LU_decompose()
{
    std::cout << "LU decomposition:\n";

    // Pętla zewnętrzna - iteracja po kolumnach macierzy (etapy eliminacji)
    // j - numer aktualnej kolumny/etapu eliminacji
    for(int j = 0; j < N; j++)
    {
        // Pivoting: find max element in column if diagonal element is zero
        // Częściowy wybór elementu podstawowego (partial pivoting):
        // Jeśli element diagonalny (pivot) jest równy zero, wyszukiwany jest
        // maksymalny co do wartości bezwzględnej element w bieżącej kolumnie
        // poniżej aktualnego wiersza, aby uniknąć dzielenia przez zero
        // i poprawić stabilność numeryczną algorytmu
        if(A[rowsIndexes[j]][j] == 0.0)
        {
            // Inicjalizacja: zakładamy, że maksimum jest w bieżącym wierszu
            int maxRow = j;
            double max = fabs(A[rowsIndexes[j]][j]);

            // Przeszukiwanie kolumny j w wierszach poniżej wiersza j
            // w celu znalezienia elementu o największej wartości bezwzględnej
            for(int i = j + 1; i < N; i++)
            {
                // Obliczenie wartości bezwzględnej elementu A[rowsIndexes[i]][j]
                double a = fabs(A[rowsIndexes[i]][j]);

                // Aktualizacja maksimum i indeksu wiersza z maksymalnym elementem
                if(a > max)
                {
                    max = a;
                    maxRow = i;
                }
            }

            // Swap row indexes to perform partial pivoting
            // Zamiana indeksów wierszy w tablicy rowsIndexes (permutacja wierszy).
            // Dzięki temu nie wykonujemy kosztownej fizycznej zamiany całych wierszy
            // w macierzy, a jedynie zamieniamy logiczne odniesienia do nich.
            // To realizuje macierz permutacji P w dekompozycji PA = LU
            int tmp = rowsIndexes[j];
            rowsIndexes[j] = rowsIndexes[maxRow];
            rowsIndexes[maxRow] = tmp;
        }

        // Gaussian elimination
        // Eliminacja Gaussa dla j-tej kolumny:
        // Dla wszystkich wierszy poniżej j-tego wiersza wykonujemy operację
        // wiersz[w] = wiersz[w] - mnożnik * wiersz[j]
        // w celu wyzerowania elementów poniżej pivota w j-tej kolumnie
        for(int w = j + 1; w < N; w++)
        {
            // Multiplier
            // Obliczenie mnożnika (współczynnika eliminacji) dla wiersza w:
            // wsp = a[w,j] / a[j,j]
            // Jest to element macierzy L poniżej przekątnej L[w,j]
            // Mnożnik określa, ile razy należy odjąć wiersz j od wiersza w,
            // aby wyzerować element w kolumnie j
            double wsp = A[rowsIndexes[w]][j] / A[rowsIndexes[j]][j];

            // Store L coefficient in A
            // Zapisanie mnożnika w miejsce elementu, który będzie wyzerowany.
            // To realizuje zapis macierzy L w dolnej części macierzy A.
            // Po zakończeniu algorytmu A[rowsIndexes[w]][j] dla w > j
            // zawiera element L[w,j] macierzy dolnotrójkątnej
            A[rowsIndexes[w]][j] = wsp;

            // Update remaining elements of the row
            // Aktualizacja elementów wiersza w na prawo od kolumny j:
            // Wykonanie operacji wiersz[w] = wiersz[w] - wsp * wiersz[j]
            // dla kolumn i > j. To jest właściwa eliminacja Gaussa.
            // Elementy w kolumnach j+1, j+2, ..., N-1 są redukowane
            // zgodnie z obliczonym mnożnikiem
            for(int i = j + 1; i < N; i++)
            {
                // Aktualizacja: a[w,i] = a[w,i] - wsp * a[j,i]
                // Tworzy to macierz górnotrójkątną U w górnej części macierzy A
                A[rowsIndexes[w]][i] -= A[rowsIndexes[j]][i] * wsp;
            }
        }

        std::cout << "Matrix A after column elimination: " << j << ":\n";
        print_matrix_A();
    }
}

//**********************************************************
//* Solves Ly = b and Ux = y using the LU-decomposed matrix.
//*
//* Procedura rozwiązująca układ równań Ax = b przy użyciu dekompozycji LU.
//* Wykorzystuje fakt, że Ax = b można zapisać jako LUx = b, co rozwiązujemy w dwóch etapach:
//* 1) Ly = b  - podstawienie w przód (forward substitution) dla układu dolnotrójkątnego
//* 2) Ux = y  - podstawienie wstecz (backward substitution) dla układu górnotrójkątnego
//* Operuje wyłącznie na wektorze b i wykorzystuje wyniki procedury LU_decompose()
//**********************************************************
void LU_solve()
{
    // Wektor pomocniczy y do przechowania rozwiązania pośredniego Ly = b
    double y[5] = {0, 0, 0, 0, 0};

    // Forward substitution: solve Ly = b
    // Podstawienie w przód - rozwiązanie układu równań Ly = b,
    // gdzie L jest macierzą dolnotrójkątną z jedynkami na przekątnej.
    // Układ rozwiązujemy od góry do dołu (od pierwszego równania do ostatniego)
    for(int i = 0; i < N; i++)
    {
        // Inicjalizacja y[i] wartością b z uwzględnieniem permutacji wierszy.
        // rowsIndexes[i] wskazuje, który oryginalny element b odpowiada
        // i-temu wierszowi po permutacji
        y[i] = b[rowsIndexes[i]];

        // Odejmowanie wpływu wcześniej obliczonych wartości y[j] (dla j < i)
        // pomnożonych przez odpowiednie elementy macierzy L.
        // Realizuje wzór: y[i] = b[i] - suma(L[i,j] * y[j]) dla j = 0..i-1
        for(int j = 0; j < i; j++)
        {
            // L[i,j] jest zapisane w A[rowsIndexes[i]][j] dla i > j
            // Odejmujemy wkład j-tej zmiennej y[j] z równania i-tego
            y[i] -= A[rowsIndexes[i]][j] * y[j];
        }
        // Po zakończeniu pętli wewnętrznej y[i] zawiera rozwiązanie i-tego
        // równania układu Ly = b (ponieważ L[i,i] = 1, nie dzielimy)
    }

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\ny after solving Ly = b\n";
    for(int i = 0; i < N; ++i)
    {
        std::cout << "y[" << rowsIndexes[i] << "] = " << y[rowsIndexes[i]] << std::endl;
    }

    // Back substitution: solve Ux = y
    // Podstawienie wstecz - rozwiązanie układu równań Ux = y,
    // gdzie U jest macierzą górnotrójkątną z elementami niezerowymi na przekątnej.
    // Układ rozwiązujemy od dołu do góry (od ostatniego równania do pierwszego)
    for(int i = 4; i >= 0; i--)
    {
        x[i] = y[i];

        // Odejmowanie wpływu już obliczonych wartości x[j] (dla j > i)
        // pomnożonych przez odpowiednie elementy macierzy U.
        // Realizuje wzór: x[i] = (y[i] - suma(U[i,j] * x[j])) / U[i,i] dla j = i+1..N-1
        for(int j = i + 1; j < N; j++)
        {
            // U[i,j] jest zapisane w A[rowsIndexes[i]][j] dla j >= i
            // Odejmujemy wkład j-tej zmiennej x[j], która została już obliczona
            // (iterujemy od końca, więc x[j] dla j > i są już znane)
            x[i] -= A[rowsIndexes[i]][j] * x[j];
        }

        // Podzielenie przez element diagonalny U[i,i] w celu otrzymania x[i].
        // Realizuje normalizację równania: x[i] = (prawa_strona) / U[i,i]
        // Po tej operacji x[i] zawiera finalną wartość i-tej niewiadomej
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