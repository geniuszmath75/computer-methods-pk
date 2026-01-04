#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

#define RESULTS 1
#define ERRORS 2

/**
 * Struktura przechowująca dane do zapisu:
 * - x        : punkt siatki
 * - u_x      : rozwiązanie analityczne U(x)
 * - uc_x     : rozwiązanie numeryczne metodą trzypunktową
 * - us_x     : rozwiązanie numeryczne metodą strzałów
 * - log10_h  : log10 kroku siatki
 * - log10_err_3pkt   : log10 maksymalnego błędu metody trzypunktowej
 * - log10_err_shoot  : log10 maksymalnego błędu metody strzałów
 */
struct ResultsRecord
{
    double x;
    double u_x;
    double uc_x;
    double us_x;
    double log10_h;
    double log10_err_3pkt;
    double log10_err_shoot;
};

/**
 * Funkcja realizuje fazę eliminacji w przód algorytmu Thomasa
 * dla układu równań z macierzą trójdiagonalną.
 *
 * Macierz A zapisana jest w postaci:
 * A[3*i]   – element poddiagonalny (a_i)
 * A[3*i+1] – element diagonalny    (b_i)
 * A[3*i+2] – element naddiagonalny (c_i)
 *
 * Metoda Thomasa jest używana, ponieważ powstający układ równań
 * po dyskretyzacji ODE jest trójdiagonalny.
 */
void factor_thomas(int n, double* A)
{
    for (int i = 1; i < n; i++)
    {
        // Modyfikacja współczynnika poddiagonalnego (a_i / b_{i-1})
        A[3*i] /= A[3*i - 2];

        // Aktualizacja współczynnika diagonalnego
        A[3*i + 1] -= A[3*i] * A[3*i - 1];
    }
}

/**
 * Funkcja wykonuje podstawianie w przód i wstecz
 * algorytmu Thomasa w celu wyznaczenia wektora rozwiązania x.
 *
 * Parametry:
 * - M : zmodyfikowana macierz po eliminacji w przód
 * - x : prawa strona układu, nadpisywana przez rozwiązanie
 */
void solve_thomas(int n, double *M, double *x)
{
    // Podstawianie w przód – modyfikacja prawej strony
    for (int i = 1; i < n; i++)
    {
        x[i] -= M[3*i] * x[i - 1];
    }

    // Wyznaczenie ostatniej niewiadomej
    x[n - 1] /= M[3 * (n - 1) + 1];

    // Podstawianie wstecz – rozwiązanie całego układu
    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = (x[i] - M[3*i + 2] * x[i + 1]) / M[3*i + 1];
    }

}

/**
 * Funkcja rozwiązuje zagadnienie brzegowe metodą trzypunktową
 * (konwencjonalna dyskretyzacja drugiego rzędu).
 *
 * Równanie:
 * U'' + 2U' - 4U + (x^3) / 2 = 0
 *
 * Warunki brzegowe:
 * U(0) = ua, U(1) = ub
 *
 * Wykorzystuje klasyczny schemat różnicowy drugiego rzędu
 * oraz algorytm Thomasa do rozwiązania powstałego układu liniowego.
 */
void solve_conventional(int n, double ua, double ub, double* xi, double* ui)
{
    // Krok siatki jednorodnej
    double h = 1.0 / (n - 1);

    // Tablica współczynników macierzy trójdiagonalnej
    double* A = new double[3*n];

    // Wektor prawej strony oraz rozwiązania
    double* u = new double[n];

    // Warunek brzegowy w x = 0
    A[0] = 0.0;
    A[1] = 1.0;
    A[2] = 0.0;
    u[0] = ua;

    // Węzły wewnętrzne – dyskretyzacja równania różniczkowego
    for(int i = 1; i < n - 1; ++i)
    {
        // Współrzędna punktu siatki
        xi[i] = i * h;

        // Współczynniki schematu różnicowego
        A[3*i] = 1.0 - h;                   // współczynnik przy u_{i-1}
        A[3*i + 1] = -2.0 - 4.0 * h * h;    // współczynnik przy u_i
        A[3*i + 2] = 1.0 + h;               // współczynnik przy u_{i+1}

        // Prawa strona wynikająca z członu x^3
        u[i] = -0.5 * h * h * xi[i] * xi[i] * xi[i];
    }

    // Warunek brzegowy w x = 1
    A[3*(n-1)] = 0.0;  
    A[3*(n-1) + 1] = 1.0;
    A[3*(n-1) + 2] = 0.0;
    u[n-1] = ub;

    // Rozwiązanie układu trójdiagonalnego
    factor_thomas(n, A);
    solve_thomas(n, A, u);

    // Przepisanie wyników do tablic wyjściowych
    for (int i = 0; i < n; ++i)
    {
        ui[i] = u[i];
        xi[i] = i * h;
    }
}

/**
 * Funkcja rozwiązuje zagadnienie brzegowe metodą strzałów.
 *
 * Zagadnienie brzegowe sprowadzane jest do zagadnienia początkowego
 * poprzez wprowadzenie niewiadomej pochodnej U'(0) = p.
 *
 * Parametr p dobierany jest iteracyjnie metodą siecznych,
 * tak aby spełniony był warunek U(1) = ub.
 */
void solve_shooting(int n, double ua, double ub, double* xi, double* ui)
{
    // Krok siatki
    double h = 1.0 / (n - 1);

    // Przybliżenia rozwiązania w kolejnych węzłach
    double* y = new double[n];
    
    // Współczynniki schematu różnicowego (jak w metodzie konwencjonalnej)
    double alfa = 1.0 - h;
    double beta = -2.0 - 4.0 * h * h;
    double gamma = 1.0 + h;

    // Funkcja obliczająca residuum warunku brzegowego w x = 1
    auto residuum = [&](double p)
    {
        // Warunki początkowe
        y[0] = ua;
        y[1] = y[0] + h * p;

        // Integracja równania różnicowego w przód
        for(int i = 1; i < n - 1; ++i)
        {
            xi[i] = i * h;
            double b = -0.5 * h*h * xi[i] * xi[i] * xi[i];
            y[i+1] = (b - alfa * y[i-1] - beta * y[i]) / gamma;
        }

        // Odchylenie od warunku brzegowego
        return y[n-1] - ub;
    };

    // Iteracyjny dobór p – metoda siecznych
    // Celem jest znalezienie takiej wartości parametru p,
    // dla której residuum(p) ≈ 0 (warunek spełnienia równania)

    // Początkowe przybliżenie p0 – różnica granic przedziału
    double p0 = (ub - ua);

    // Drugie przybliżenie p1 – niewielkie przesunięcie względem p0
    double p1 = p0 + 1.0;


    // Wartości funkcji (residuum) w punktach początkowych
    double f0 = residuum(p0);
    double f1 = residuum(p1);

    // Kolejne przybliżenie rozwiązania
    double p2 = p1;

    // Główna pętla iteracyjna metody siecznych
    for(int iter = 0; iter < 100; ++iter)
    {
        // Sprawdzenie, czy mianownik we wzorze metody siecznych
        // nie jest zbyt mały (f1 ≈ f0).
        // Zbyt mała różnica prowadziłaby do niestabilności numerycznej
        // lub dzielenia przez liczbę bliską zeru.
        if (fabs(f1 - f0) < 1e-12) break;

        // Wzór metody siecznych:
        p2 = p1 - f1 * (p1 - p0) / (f1 - f0);

        // Obliczenie residuum dla nowego przybliżenia
        double f2 = residuum(p2);
        
        // Warunek zakończenia iteracji:
        // 1) |p2 - p1| < 1e-8  → kolejne przybliżenia są bardzo blisko siebie
        //    (stabilizacja rozwiązania)
        // 2) |f2| < 1e-8       → residuum bliskie zeru
        //    (równanie spełnione z zadaną dokładnością)
        if (fabs(p2 - p1) < 1e-12 && fabs(f2) < 1e-12) break;
        
        p0 = p1;
        f0 = f1;
        p1 = p2;
        f1 = f2;
    }

    // Przepisanie rozwiązania końcowego
    for(int i = 0; i < n; ++i) {
        xi[i] = i * h;
        ui[i] = y[i];
    }
}

/**
 * Funkcja zwraca rozwiązanie analityczne U(x)
 * podane w treści zadania.
 *
 * Wykorzystywana jako rozwiązanie referencyjne
 * do obliczania błędów bezwzględnych metod numerycznych.
 */
double U_analytic(double x)
{
    const double sqrt5 = sqrt(5.0);
    double term2 = 95.0 * exp((-1.0 - sqrt5) * (-1.0 + x));
    double term3 = 55.0 * exp((-1.0 + sqrt5) * x);
    double term4 = 95.0 * exp( (1.0 + sqrt5) + (-1.0 + sqrt5) * x);
    double term5 = -55.0 * exp(2.0 * sqrt5 - (1.0 + sqrt5) * x);
    double poly1 = 2.0 * x * (6.0 + x * (3.0 + 2.0 * x));
    double term6 = exp(2.0 * sqrt5) * (9.0 + 2.0 * x * (6.0 + x*(3.0 + 2.0*x)));
    double bracket = 9.0 - term2 + term3 + term4 + term5 + poly1 - term6;
    double coth5 = cosh(sqrt5) / sinh(sqrt5);
    return -((bracket * (-1.0 + coth5)) / 64.0);
}

/**
 * Funkcja zapisująca dane do pliku
 */
void save_results(string fileName, const vector<ResultsRecord>& dataRecords, int typeId )
{
    ofstream fout(fileName);
    if(!fout.is_open())
    {
        cerr << "[ERROR]: Cannot open output file: " << fileName << "\n";
        return;
    }

    size_t i = 0;

    if (typeId == ERRORS)
    {
        fout << "#log10(h)       log10(err_3pkt)      log10(err_strzal)\n";
    }
    if (typeId == RESULTS)
    {
        fout << "#x         U(x)         UC(x)       US(x)\n";
    }
    fout << fixed << setprecision(12);

    while (i < dataRecords.size())
    {
        if (typeId == ERRORS)
        {
            fout << dataRecords[i].log10_h << "\t\t"
             << dataRecords[i].log10_err_3pkt << "\t\t\t"
             << dataRecords[i].log10_err_shoot << endl;
        }
        if(typeId == RESULTS)
        {
            fout << dataRecords[i].x << "\t"
                 << dataRecords[i].u_x << "\t"
                 << dataRecords[i].uc_x << "\t" 
                 << dataRecords[i].us_x << endl;
        }
        
        i++;
    }

    fout.close();
    cout << "Results saved to file: " << fileName << endl;
}


/**
 * Funkcja oblicza nachylenie prostej regresji liniowej
 * y = a x + b metodą najmniejszych kwadratów.
 *
 * Nachylenie a interpretowane jest jako
 * doświadczalny rząd dokładności metody.
 */
double linear_regression_slope(
    const vector<double>& x,
    const vector<double>& y)
{
    int n = x.size();
    // Suma x, y, x^2, y^2
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;

    for (int i = 0; i < n; ++i)
    {
        sx  += x[i];
        sy  += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
    }

    // (y2 - y1) / (x2 - x1)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}


/**
 * Funkcja główna realizuje kolejne etapy zadania:
 * 1. Obliczenie rozwiązań numerycznych
 * 2. Porównanie z rozwiązaniem analitycznym
 * 3. Analizę błędów i rzędu dokładności
 */
int main()
{
    // Wektory do przechowywania wyników i błędów
    vector<ResultsRecord> resultRecords;
    vector<ResultsRecord> errorRecords;

    // Liczba węzłów do wypisania przykładowych wyników
    int N = 10;

    // Tablice punktów siatki oraz rozwiązań
    double* xi_1 = new double[N];
    double* uc_1 = new double[N];
    double* us_1 = new double[N];

    // Rozwiązanie zadania metodą trzypunktową i strzałów
    solve_conventional(N, 2.0, -2.0, xi_1, uc_1);
    solve_shooting(N, 2.0, -2.0, xi_1, us_1);

    cout << "x        U(x)        UC(x)       US(x)\n";
    cout << fixed << setprecision(8);

    // Wypisanie porównania trzech rozwiązań
    for (int i = 0; i < N; i++)
    {   
        double ua = U_analytic(xi_1[i]);
        ResultsRecord resultRecord;
        resultRecord.x = xi_1[i];
        resultRecord.u_x = ua;
        resultRecord.uc_x = uc_1[i];
        resultRecord.us_x = us_1[i];
        resultRecords.push_back(resultRecord);
        cout << setw(8) << xi_1[i] << "  " << ua << "  " << uc_1[i] << "  " << us_1[i] << "\n";
    }

     // Zapis wyników punktowych
    save_results("wyniki_x_ux.txt", resultRecords, RESULTS);

    // Analiza błędu dla różnych kroków siatki
    for (int n = 10; n <= 1000000; n *= 1.25)
    {
        double h = 1.0 / (n - 1);
        double* xi = new double[n];
        double* uc = new double[n];
        double* us = new double[n];

        solve_conventional(n, 2.0, -2.0, xi, uc);
        solve_shooting(n, 2.0, -2.0, xi, us);

        double max_err_c = 0.0;
        double max_err_s = 0.0;

        for (int i = 0; i < n; i++)
        {
            double err_c = fabs(uc[i] - U_analytic(xi[i]));
            double err_s = fabs(us[i] - U_analytic(xi[i]));

            if (err_c > max_err_c)
                max_err_c = err_c;
            if (err_s > max_err_s)
                max_err_s = err_s;
        }

        ResultsRecord errorsRecord;
        errorsRecord.log10_h = log10(h);
        errorsRecord.log10_err_3pkt = log10(max_err_c);
        errorsRecord.log10_err_shoot = log10(max_err_s);
        errorRecords.push_back(errorsRecord);
    }

    // Zapis danych do analizy rzędu dokładności
    save_results("wyniki_h_errors.txt", errorRecords, ERRORS);

    // =========================================================
    // ANALIZA RZĘDU DOKŁADNOŚCI I BŁĘDÓW MASZYNOWYCH
    // =========================================================

    vector<double> logh, loge_c, loge_s;

    for (const auto& rec : errorRecords)
    {
        // Pomijamy punkty z -inf (zbyt mały błąd) i zbyt małym krokiem h (wpływ błędów maszynowych)
        if (isfinite(rec.log10_err_3pkt) && isfinite(rec.log10_err_shoot) && rec.log10_h >= -4.2)
        {
            logh.push_back(rec.log10_h);
            loge_c.push_back(rec.log10_err_3pkt);
            loge_s.push_back(rec.log10_err_shoot);
        }
    }

    double order_conv = linear_regression_slope(logh, loge_c);
    double order_shot = linear_regression_slope(logh, loge_s);

    cout << "\n=============================================\n";
    cout << "DOŚWIADCZALNE RZĘDY DOKŁADNOŚCI\n";
    cout << "=============================================\n";
    cout << "Metoda trzypunktowa:  p ≈ " << order_conv << " (teoria: 2)\n";
    cout << "Metoda strzałów:      p ≈ " << order_shot << " (teoria: 2)\n";

    // ---------------------------------------------------------
    // Identyfikacja wpływu błędów maszynowych
    // ---------------------------------------------------------
    auto min_c = min_element(loge_c.begin(), loge_c.end());
    auto min_s = min_element(loge_s.begin(), loge_s.end());

    int idx_c = distance(loge_c.begin(), min_c);
    int idx_s = distance(loge_s.begin(), min_s);

    cout << "\n=============================================\n";
    cout << "WPŁYW BŁĘDÓW MASZYNOWYCH\n";
    cout << "=============================================\n";
    cout << "Metoda trzypunktowa:\n";
    cout << "  minimalny błąd przy h ≈ "
         << logh[idx_c] << "\n";

    cout << "Metoda strzałów:\n";
    cout << "  minimalny błąd przy h ≈ "
         << logh[idx_s] << "\n";


    return 0;
}