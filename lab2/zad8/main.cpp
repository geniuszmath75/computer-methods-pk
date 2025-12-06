#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
#include <fstream>

/*
 Funkcja safe_log10:
 Zwraca log10(e), ale zabezpiecza przed próbą obliczenia logarytmu
 z liczby 0 lub zbyt małej (poniżej najmniejszej normalnej liczby
 danego typu). W takim przypadku zwracamy NaN, aby nie psuć wykresów.
*/
template <typename T>
T safe_log10(T e)
{
    // najmniejsza normalna dodatnia liczba w danym typie
    const T tiny = std::numeric_limits<T>::min();
    if (e <= tiny || !std::isfinite(e))
        return std::numeric_limits<T>::quiet_NaN();
    else
        return std::log10(e);
}

/*
 Stała PI w długiej precyzji (long double)
*/
namespace math_const
{
    constexpr long double PI = 3.141592653589793238462643383279502884L;
}


/* Funkcja f(x) = cos(x) */
template <typename T>
inline T f(T x)
{
    return std::cos(x);
}

/* Dokładna pochodna f(x): f'(x) = -sin(x) */
template <typename T>
inline T df_exact(T x)
{
    return -std::sin(x);
}

/* Dwupunktowa pochodna w przód */
template <typename T>
inline T fd_forward2(T (*func)(T), T x, T h)
{
    return (func(x + h) - func(x)) / h;
}

/* Trzypunktowa pochodna w przód */
template <typename T>
inline T fd_forward3(T (*func)(T), T x, T h)
{
    return (-T(3) / T(2) * func(x) + T(2) * func(x + h) - T(1) / T(2) * func(x + T(2) * h)) / h;
}

/* Dwupunktowa pochodna w tył */
template <typename T>
inline T fd_backward2(T (*func)(T), T x, T h)
{
    return (func(x) - func(x - h)) / h;
}

/* Trzypunktowa pochodna w tył */
template <typename T>
inline T fd_backward3(T (*func)(T), T x, T h)
{
    return (T(1) / T(2) * func(x - T(2) * h) - T(2) * func(x - h) + T(3) / T(2) * func(x)) / h;
}

/* Dwupunktowa pochodna centralna */
template <typename T>
inline T fd_central2(T (*func)(T), T x, T h)
{
    return (func(x + h) - func(x - h)) / (T(2) * h);
}

/*
 Funkcja zapisująca dane do pliku.
 Format pliku:
 log10(h) log10(e_f2) log10(e_f3) log10(e_c2) log10(e_b2) log10(e_b3)
*/
template <typename T>
void save_to_file(
    const std::string &filename,
    const std::vector<T> &logh,
    const std::vector<T> &ef2,
    const std::vector<T> &ef3,
    const std::vector<T> &ec2,
    const std::vector<T> &eb2,
    const std::vector<T> &eb3)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Błąd: nie można otworzyć pliku " << filename << " do zapisu.\n";
        return;
    }

    file << "# log10(h) log10(e_f2) log10(e_f3) log10(e_c2) log10(e_b2) log10(e_b3)\n";

    for (size_t i = 0; i < logh.size(); ++i)
    {
        file << std::setprecision(10)
             << logh[i] << " "
             << ef2[i] << " "
             << ef3[i] << " "
             << ec2[i] << " "
             << eb2[i] << " "
             << eb3[i] << "\n";
    }

    std::cout << "Zapisano wyniki do pliku: " << filename << "\n";
}

/*
 Funkcja wykonuje eksperyment numeryczny dla danego typu T:
 - oblicza pochodne w trzech punktach: x0, xc, xn
 - liczy błędy dla różnych metod różnicowych
 - wypisuje wyniki
 - zapisuje dane do pliku

 label – label typu (np. "double")
*/
template <typename T>
void run_for_type(const std::string &label)
{
    const T PI_2 = static_cast<T>(math_const::PI / 2.0L);
    const T x0 = T(0);          // początek przedziału
    const T xn = PI_2;          // koniec przedziału
    const T xc = PI_2 / T(2);   // środek przedziału

    std::cout << "# type = " << label << "\n";

    // wektory do zapisu do pliku
    std::vector<T> logh, ef2, ef3, ec2, eb2, eb3;

    std::cout << std::setw(13) << "h"
              << std::setw(14) << "|e_f2(x0)|"
              << std::setw(14) << "|e_f3(x0)|"
              << std::setw(14) << "|e_c2(xc)|"
              << std::setw(14) << "|e_b2(xn)|"
              << std::setw(14) << "|e_b3(xn)|"
              << "\n";

    for (int k = 0; k <= 50; ++k)
    {
        // krok sieci: h = (pi/2) * 2^{-k}
        T h = std::ldexp(T(PI_2), -k);

        if (T(2) * h > PI_2)
        {
            continue;
        }

        const T de_x0 = df_exact(x0);
        const T de_xc = df_exact(xc);
        const T de_xn = df_exact(xn);

        // obliczanie przybliżeń różnicowych
        const T d_f2_x0 = fd_forward2(f<T>, x0, h);
        const T d_f3_x0 = fd_forward3(f<T>, x0, h);
        const T d_c2_xc = fd_central2(f<T>, xc, h);
        const T d_b2_xn = fd_backward2(f<T>, xn, h);
        const T d_b3_xn = fd_backward3(f<T>, xn, h);

        // błędy bezwzględne
        const T e_f2_x0 = std::abs(d_f2_x0 - de_x0);
        const T e_f3_x0 = std::abs(d_f3_x0 - de_x0);
        const T e_c2_xc = std::abs(d_c2_xc - de_xc);
        const T e_b2_xn = std::abs(d_b2_xn - de_xn);
        const T e_b3_xn = std::abs(d_b3_xn - de_xn);

        // zapisujemy log10 błędów
        logh.push_back(std::log10(h));
        ef2.push_back(safe_log10(e_f2_x0));
        ef3.push_back(safe_log10(e_f3_x0));
        ec2.push_back(safe_log10(e_c2_xc));
        eb2.push_back(safe_log10(e_b2_xn));
        eb3.push_back(safe_log10(e_b3_xn));

        std::cout << std::scientific << std::setprecision(6);

        std::cout << std::log10(static_cast<T>(h)) << " "
                  << safe_log10(static_cast<T>(e_f2_x0)) << " "
                  << safe_log10(static_cast<T>(e_f3_x0)) << " "
                  << safe_log10(static_cast<T>(e_c2_xc)) << " "
                  << safe_log10(static_cast<T>(e_b2_xn)) << " "
                  << safe_log10(static_cast<T>(e_b3_xn)) << "\n";
    }
    
    // zapisz wyniki do pliku CSV
    save_to_file(
        "wyniki_" + label + ".txt",
        logh, ef2, ef3, ec2, eb2, eb3
    );
}

int main()
{
    run_for_type<double>("double");
    run_for_type<long double>("long_double");
    return 0;
}