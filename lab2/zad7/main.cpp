#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Stała określająca rozmiar układu równań (5 równań, 5 niewiadomych)
const int N = 5;
const int ITERATIONS = 50; // Max. liczba iteracji


// Wygodne aliasy na wektor i macierz
typedef vector<double> Vec;
typedef vector<Vec> Matrix;

// ============================================================================
// Funkcja obliczająca residuum r = b - A*x
// 
// ============================================================================
Vec compute_residual(const Matrix &A, const Vec &x, const Vec &b)
{
    Vec r(N);
    for (int i = 0; i < N; ++i)     // przejście po kolejnych równaniach
    {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) // liczymy sumę A[i][j] * x[j]
        {
            sum += A[i][j] * x[j];
        }
        r[i] = b[i] - sum;          // residuum dla równania i-tego
    }
    return r;
}

// ============================================================================
// Funkcja wypisująca wyniki pojedynczej iteracji
// - x: aktualne przybliżenie
// - err: estymator błędu |x_new - x_old|
// - res: residuum r = b - A*x
// 
// ============================================================================
void print_iteration(int iteration, const Vec &x, const Vec &err, const Vec &res)
{
    cout << "Iteration " << setw(3) << fixed << setprecision(6) << iteration << ":\n" << "x = [";
    for (int i = 0; i < N; ++i)
    {
        cout << setw(10) << x[i] << (i + 1 < N ? ", " : " ],\n");
    }
    cout << "err = [" << setprecision(6);
    for (int i = 0; i < N; ++i)
    {
        cout << setw(10) << err[i] << (i + 1 < N ? ", " : " ],\n");
    }
    cout << "res = [";
    for (int i = 0; i < N; ++i)
    {
        cout << setw(10) << res[i] << (i + 1 < N ? ", " : " ],\n");
    }
    cout << '\n';
}

// ============================================================================
// Funkcja sprawdzająca kryterium zbieżności
// - sprawdza dwa warunki:
//   1) przyrost iteracji:   err <= tol_x
//   2) residuum układu:     |r| <= tol_r
//
// ============================================================================
bool is_converged(const Vec &err, const Vec &tol_x, const Vec &res, const Vec &tol_r)
{
    for (int i = 0; i < N; ++i)
    {
        // jeśli którykolwiek warunek nie jest spełniony — brak zbieżności
        if (err[i] > tol_x[i] || fabs(res[i]) > tol_r[i])
            return false;
    }
    return true;
}

// ============================================================================
// METODA JACOBIEGO
// x_new obliczane wyłącznie z x_old (nie używa zaktualizowanych wartości)
// 
// ============================================================================
void jacobi_method(const Matrix &A, const Vec &b, Vec x0, const Vec &tol_x, const Vec &tol_r)
{
    cout << "=== JACOBI METHOD ===\n";
    Vec x_old = x0;     // poprzednie przybliżenie
    Vec x_new(N);       // nowe przybliżenie
    Vec err(N), res(N); // estymator błędu i residuum

    int iteration = 0;

    do
    {
        ++iteration;
        // --------------------------------------------------------------------
        // Iteracja Jacobiego
        // --------------------------------------------------------------------
        for (int i = 0; i < N; ++i)
        {
            double sigma = 0.0;

            // oblicz sumę po wszystkich j != i
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    sigma += A[i][j] * x_old[j];
                }
            }

            // wzór Jacobiego
            x_new[i] = (b[i] - sigma) / A[i][i];

            // błąd = różnica między przybliżeniami
            err[i] = fabs(x_new[i] - x_old[i]);
        }

        // obliczamy residuum r = b - A*x_new
        res = compute_residual(A, x_new, b);

        // wypisujemy bieżące wyniki
        print_iteration(iteration, x_new, err, res);

        // przygotowanie do następnej iteracji: x_old = x_new
        x_old = x_new;
    } while (!is_converged(err, tol_x, res, tol_r) && iteration < ITERATIONS);    // sprawdzamy warunek stopu
    cout << "Iteration number: " << iteration << "\n\n";
}

// ============================================================================
// METODA GAUSSA-SEIDELA
// Używa zaktualizowanych wartości x[n] natychmiast w kolejnych obliczeniach
// 
// ============================================================================
void gauss_seidel_method(const Matrix &A, const Vec &b, Vec x0, const Vec &tol_x, const Vec &tol_r)
{
    cout << "=== GAUSS-SEIDEL METHOD ===\n";
    Vec x_old = x0;     // poprzednie przybliżenie (tu pomocniczo)
    Vec x_new = x0;     // aktualne (aktualizowane w miejscu)
    Vec err(N), res(N);

    int iteration = 0;

    do
    {
        ++iteration;

        // --------------------------------------------------------------------
        // Iteracja Gaussa-Seidela
        // --------------------------------------------------------------------
        for (int i = 0; i < N; ++i)
        {
            double sigma = 0.0;

            // sumujemy po wszystkich j != i
            // UWAGA: x_new zawiera już część wartości z bieżącej iteracji!
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                {
                    sigma += A[i][j] * x_new[j];
                }
            }

            double xi = (b[i] - sigma) / A[i][i];
            
            // obliczamy estymator błędu
            err[i] = fabs(xi - x_new[i]);

            // zapisujemy nową wartość
            x_new[i] = xi;
        }
        // obliczamy residuum
        res = compute_residual(A, x_new, b);
        
        print_iteration(iteration, x_new, err, res);
    } while (!is_converged(err, tol_x, res, tol_r) && iteration < ITERATIONS);
    cout << "Iteration number: " << iteration << "\n\n";
}


// ============================================================================
// METODA SOR (Successive Over-Relaxation)
// Wersja GS ze współczynnikiem relaksacji ω
// omega = 1 -> metoda Gaussa-Seidela
// 0 < ω < 1 -> niedorelaksacja (wolniej, stabilniej)
// 
// ============================================================================
void sor_method(const Matrix &A, const Vec &b, Vec x0, double omega, const Vec &tol_x, const Vec &tol_r)
{
    cout << "=== SOR METHOD (omega = " << omega << ") ===\n";
    Vec x_new = x0;
    Vec err(N), res(N);

    int iteration = 0;

    do
    {
        ++iteration;

        Vec x_old = x_new;  // zachowujemy poprzednie wartości (do liczenia błędu)
        
        // --------------------------------------------------------------------
        // Iteracja SOR
        // --------------------------------------------------------------------
        for (int i = 0; i < N; ++i)
        {
            double sigma = 0.0;

            // suma po pozostałych elementach
            for (int j = 0; j < N; ++j)
            {
                if (j != i)
                {
                    sigma += A[i][j] * x_new[j];
                }
            }

            // wzór SOR:
            double xi = (1 - omega) * x_old[i] + (omega / A[i][i]) * (b[i] - sigma);

            // estymator błędu
            err[i] = fabs(xi - x_old[i]);
            
            x_new[i] = xi;
        }

        // residuum
        res = compute_residual(A, x_new, b);
        
        print_iteration(iteration, x_new, err, res);

    } while (!is_converged(err, tol_x, res, tol_r) && iteration < ITERATIONS);
    cout << "Iteration number: " << iteration << "\n\n";
}

int main()
{
    // Definicja macierzy A 
    Matrix A = {
        { 50,  5,  4,   3,  2 },
        {  1, 40,  1,   2,  3 },
        {  4,  5, 30,  -5, -4 },
        { -3, -2, -1, -20,  0 },
        {  1,  2,  3,   4, 30 }
    };

    // Wektor b
    Vec b = { 140, 67, 62, 89, 153 };
    
    // Przybliżenie początkowe
    Vec x0 = { 6, 6, 6, 6, 6 };

    // Tolerancje błędu i residuum (np. 1e-6 dla każdego elementu)
    Vec tol_x = { 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
    Vec tol_r = { 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
    double omega = 0.5;

    jacobi_method(A, b, x0, tol_x, tol_r);
    gauss_seidel_method(A, b, x0, tol_x, tol_r);
    sor_method(A, b, x0, omega, tol_x, tol_r);

    return 0;
}