#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

/*
 Struktura przechowująca rozwiązania równania różniczkowego
 w kolejnych punktach siatki czasowej.
 Każdy rekord odpowiada jednemu czasowi t_k.
*/
struct ResultRecord
{
    double t;
    double y_exact;
    double y_bme;
    double y_pme;
    double y_trap;
};

/*
 Struktura przechowująca dane do analizy błędów globalnych.
 Wartości zapisane są w skali logarytmicznej, co odpowiada
 wykresowi log10|błąd| w funkcji log10(dt).
*/
struct ErrorRecord
{
    double log10_dt;
    double log10_err_bme;
    double log10_err_pme;
    double log10_err_trapez;
};


/*
 Funkcja zwracająca współczynnik a(t) występujący w równaniu różniczkowym:
 y'(t) + a(t)(y(t) - 1) = 0.
*/
double a_coeff(double t)
{
    return (100.0 * t + 10.0) / (t + 1.0);
}

/*
 Funkcja zwracająca dokładne (analityczne) rozwiązanie równania różniczkowego.
*/
double y_exact(double t)
{
    return 1.0 + pow(1.0 + t, 90.0) * exp(-100.0 * t);
}

/*
 Funkcja realizująca jeden krok czasowy trzema metodami numerycznymi:
 - Eulerem bezpośrednim,
 - Eulerem pośrednim,
 - metodą trapezów.
 
 Argumenty:
 t        – aktualny czas t_k,
 dt       – krok czasowy Δt,
 y_*_k    – wartości rozwiązania w chwili t_k,
 y_*_k1   – wartości rozwiązania w chwili t_{k+1}.
*/
void step_methods(double t, double dt,
                  double y_bme_k, double y_pme_k, double y_trap_k,
                  double& y_bme_k1, double& y_pme_k1, double& y_trap_k1)
{
    // Współczynnik a(t_k) oraz a(t_{k+1})
    double ak = a_coeff(t);
    double ak1 = a_coeff(t + dt);

    // 1) Euler bezpośredni (BME)
    // y_{k+1} = y_k + Δt * f(t_k, y_k)
    // gdzie f(t,y) = -a(t)(y-1)
    y_bme_k1 = y_bme_k + dt * (-ak * (y_bme_k - 1.0));

    // 2) Euler pośredni (PME)
    // y_{k+1} = y_k + Δt * f(t_{k+1}, y_{k+1})
    // po przekształceniu otrzymujemy wzór jawny
    y_pme_k1 = ( y_pme_k + dt * ak1 ) / ( 1.0 + dt * ak1 );

    // 3) Metoda trapezów
    // y_{k+1} = y_k + Δt/2 * ( f(t_k,y_k) + f(t_{k+1},y_{k+1}) )
    // równanie rozwiązane względem y_{k+1}
    y_trap_k1 = ( y_trap_k - (dt / 2.0) * ak * (y_trap_k - 1.0) + (dt / 2.0) * ak1 ) / ( 1.0 + (dt / 2.0) * ak1 );
}


/*
 Funkcja obliczająca maksymalne błędy bezwzględne
 rozwiązań numerycznych dla danego kroku czasowego Δt.
 
 Błąd globalny definiowany jest jako:
 max_k | y_num(t_k) - y_exact(t_k) |.
*/
void compute_max_errors(
    double T, int N,
    double &err_bme, double &err_pme, double &err_trapez)
{
    double dt = T / N;

    // Warunek początkowy y(0)=2 dla wszystkich metod
    double y_fe = 2.0, y_be = 2.0, y_tr = 2.0;
    err_bme = err_pme = err_trapez = 0.0;

    // Iteracja po całym przedziale czasu
    for (int k = 0; k < N; ++k)
    {
        double t = k * dt;
        double y_fe_next, y_be_next, y_tr_next;

        // Jeden krok metod numerycznych
        step_methods(t, dt, y_fe, y_be, y_tr,
                     y_fe_next, y_be_next, y_tr_next);

        double t_next = t + dt;
        double y_ex = y_exact(t_next);

        // Aktualizacja maksymalnego błędu bezwzględnego
        err_bme = max(err_bme, fabs(y_fe_next - y_ex));
        err_pme = max(err_pme, fabs(y_be_next - y_ex));
        err_trapez = max(err_trapez, fabs(y_tr_next - y_ex));

        y_fe = y_fe_next;
        y_be = y_be_next;
        y_tr = y_tr_next;
    }
}

/*
 Funkcja zapisująca rozwiązania numeryczne i analityczne do pliku.
*/
void save_results(const string& fileName, vector<ResultRecord> &data)
{
    ofstream fout(fileName);
    if (!fout.is_open())
    {
        cerr << "[ERROR]: Cannot open output file: " << fileName << 
        endl;
        return;
    }

    fout << scientific << setprecision(8) << setw(10);
    fout << "#t\t\ty_exact\t\ty_BME\t\ty_PME\t\ty_TRAP" << endl;
    for (const auto &r : data)
    {
        fout << r.t << "\t"
             << r.y_exact << "\t"
             << r.y_bme << "\t"
             << r.y_pme << "\t"
             << r.y_trap << endl;
    }
}

/*
 Funkcja zapisująca dane błędów do osobnego pliku.
*/
void save_results(const std::string &fileName,
                  const std::vector<ErrorRecord> &data)
{
    std::ofstream fout(fileName);
    if (!fout.is_open())
    {
        std::cerr << "[ERROR] Cannot open file: " << fileName << "\n";
        return;
    }

    fout << scientific << setprecision(8) << setw(10);
    fout << "#log10(dt)\t\tlog10(err_BME)\t\tlog10(err_PME)\t\tlog10(err_TRAPEZ)" << endl;
    for (const auto &e : data)
    {
        fout << e.log10_dt << " "
             << e.log10_err_bme << " "
             << e.log10_err_pme << " "
             << e.log10_err_trapez << "\n";
    }
}

/*
 Funkcja obliczająca nachylenie prostej metodą najmniejszych kwadratów.
 W kontekście zadania nachylenie to jest doświadczalnym rzędem dokładności metody.
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


int main()
{
    double T = 5.0; // koniec przedziału czasu
    int N = 500;    // liczba kroków czasowych

    // cout << "Podaj koncowy czas T: ";
    // cin >> T;

    // cout << "Podaj liczbę kroków N:";
    // cin >> N;

    double dt = T / N;
    vector<ResultRecord> results(N + 1);
    vector<ErrorRecord> errors;

    // Warunek początkowy y(0)=2
    results[0] = { 0.0, y_exact(0.0), 2.0, 2.0, 2.0 };

    for (int k = 0; k < N; ++k)
    {
        double t = k * dt;

        step_methods(t, dt,
                    results[k].y_bme,
                    results[k].y_pme,
                    results[k].y_trap,
                    results[k + 1].y_bme,
                    results[k + 1].y_pme,
                    results[k + 1].y_trap
                );

        results[k + 1].t = t + dt;
        results[k + 1].y_exact = y_exact(t + dt);
    }

    save_results("solution.txt", results);
    cout << "Zapisano wyniki do pliku solutions.txt" << endl;

    // Lista różnych zagęszczeń siatki czasowej
    vector<int> N_values = { 10, 20, 40, 80, 155, 160, 240, 320, 480, 640, 960, 1280, 1920, 2560 };

    for (int n : N_values)
    {
        double err_bme, err_pme, err_trapez;
        compute_max_errors(T, n, err_bme, err_pme, err_trapez);

        double dt = T / n;

        errors.push_back({
            std::log10(dt),
            std::log10(err_bme),
            std::log10(err_pme),
            std::log10(err_trapez)
        });
    }

    save_results("errors.txt", errors);
    cout << "Zapisano błędy do pliku errors.txt\n";

    // ---- WYZNACZANIE RZEDOW DOKLADNOSCI ----
    vector<double> log_dt;
    vector<double> log_err_bme, log_err_pme, log_err_trapez;

    for (const auto &e : errors)
    {
        // Pomijamy wartości dt zbyt wysokie
        if(e.log10_dt < -1.5)
        {
            log_dt.push_back(e.log10_dt);
            log_err_bme.push_back(e.log10_err_bme);
            log_err_pme.push_back(e.log10_err_pme);
            log_err_trapez.push_back(e.log10_err_trapez);
        }
    }

    double p_bme = linear_regression_slope(log_dt, log_err_bme);
    double p_pme = linear_regression_slope(log_dt, log_err_pme);
    double p_trapez = linear_regression_slope(log_dt, log_err_trapez);

    std::cout << "\nDoswiadczalne rzedy dokladnosci:\n";
    std::cout << "Euler bezpośredni      p ≈ " << p_bme << "\n";
    std::cout << "Euler pośredni         p ≈ " << p_pme << "\n";
    std::cout << "Metoda trapezów        p ≈ " << p_trapez << "\n";

    return 0;
}