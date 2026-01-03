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

int nT;
void factor_thomas(int n, double* A)
{
    for (int i = 1; i < n; i++)
    {
        A[3*i] /= A[3*i - 2];
        A[3*i + 1] -= A[3*i] * A[3*i - 1];
    }
}

void solve_thomas(int n, double *M, double *x)
{
    for (int i = 1; i < n; i++)
    {
        x[i] -= M[3*i] * x[i - 1];
    }

    x[n - 1] /= M[3 * (n - 1) + 1];

    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = (x[i] - M[3*i + 2] * x[i + 1]) / M[3*i + 1];
    }

}

/**
 * 1) Trzypunktowa dyskretyzacja konwencjonalna
 */
void solve_conventional(int n, double ua, double ub, double* xi, double* ui)
{
    double h = 1.0 / (n - 1);
    double* A = new double[3*n];
    double* u = new double[n];

    A[0] = 0.0;
    A[1] = 1.0;
    A[2] = 0.0;
    u[0] = ua;

    for(int i = 1; i < n - 1; ++i)
    {
        xi[i] = i * h;
        A[3*i] = 1.0 - h;
        A[3*i + 1] = -2.0 - 4.0 * h * h;
        A[3*i + 2] = 1.0 + h;
        u[i] = -0.5 * h * h * xi[i] * xi[i] * xi[i];
    }

    A[3*(n-1)] = 0.0;  
    A[3*(n-1) + 1] = 1.0;
    A[3*(n-1) + 2] = 0.0;
    u[n-1] = ub;

    factor_thomas(n, A);
    solve_thomas(n, A, u);

    for (int i = 0; i < n; ++i)
    {
        ui[i] = u[i];
        xi[i] = i * h;
    }
}

/**
 * 2) Metoda strzałów
 */
void solve_shooting(int n, double ua, double ub, double* xi, double* ui) {
    double h = 1.0 / (n - 1);
    double* y = new double[n];
    
    double alfa = 1.0 - h;
    double beta = -2.0 - 4.0 * h * h;
    double gamma = 1.0 + h;

    auto residuum = [&](double p)
    {
        y[0] = ua;
        y[1] = y[0] + h * p;

        for(int i = 1; i < n - 1; ++i)
        {
            xi[i] = i * h;
            double b = -0.5 * h*h * xi[i] * xi[i] * xi[i];
            y[i+1] = (b - alfa * y[i-1] - beta * y[i]) / gamma;
        }

        return y[n-1] - ub;
    };

    // metoda siecznych
    double p0 = (ub - ua);
    double p1 = p0 + 1.0;
    double f0 = residuum(p0);
    double f1 = residuum(p1);
    double p2 = p1;

    for(int iter = 0; iter < 100; ++iter)
    {
        if (fabs(f1 - f0) < 1e-8) break;

        p2 = p1 - f1 * (p1 - p0) / (f1 - f0);
        double f2 = residuum(p2);
        
        if (fabs(p2 - p1) < 1e-8 && fabs(f2) < 1e-8) break;
        
        p0 = p1;
        f0 = f1;
        p1 = p2;
        f1 = f2;
    }

    for(int i = 0; i < n; ++i) {
        xi[i] = i * h;
        ui[i] = y[i];
    }
}

/**
 * Rozwiązanie analityczne U(x)
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
    fout << fixed << setprecision(8);

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
 * Wyznaczanie nachylenia regresji liniowej y = a x + b
 * (a = rząd dokładności)
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
    vector<ResultsRecord> resultRecords;
    vector<ResultsRecord> errorRecords;

    int N = 10;

    double* xi_1 = new double[N];
    double* uc_1 = new double[N];
    double* us_1 = new double[N];

    solve_conventional(N, 2.0, -2.0, xi_1, uc_1);
    solve_shooting(N, 2.0, -2.0, xi_1, us_1);

    cout << "x        U(x)        UC(x)       US(x)\n";
    cout << fixed << setprecision(8);

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

    save_results("wyniki_x_ux.txt", resultRecords, RESULTS);

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