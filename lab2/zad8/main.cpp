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
 Funkcja obliczająca lokalny doświadczalny rząd dokładności
 przybliżenia różnicowego.

 Definicja:
 p ≈ log(e₂ / e₁) / log(h₂ / h₁)

 gdzie:
 e₁, e₂ – błędy bezwzględne dla kolejnych kroków siatki,
 h₁, h₂ – odpowiadające im kroki siatki.

 W obszarze dominacji błędu obcięcia wartość p powinna
 dążyć do rzędu teoretycznego metody.
*/
template<typename T>
T experimental_order(T e1, T e2, T h1, T h2)
{
    return std::log10(e2 / e1) / std::log10(h2 / h1);
}

/*
 Funkcja obliczająca uśredniony doświadczalny rząd dokładności
 na podstawie wektora lokalnych rzędów p_i.

 Uśrednianie wykonywane jest na wybranym przedziale indeksów,
 odpowiadającym liniowemu fragmentowi wykresu log|e(h)| vs log(h),
 w którym nie występuje jeszcze wpływ błędów maszynowych.

 Wartości NaN i Inf są pomijane.
*/
template<typename T>
T average_experimental_order(
    const std::vector<T>& orders,
    size_t begin,
    size_t end
)
{
    T sum = static_cast<T>(0);
    size_t count = 0;

    for (size_t i = begin; i < end && i < orders.size(); ++i)
    {
        if(std::isfinite(orders[i]))
        {
            sum += orders[i];
            ++count;
        }
    }

    if(count == 0)
    {
        return std::numeric_limits<T>::quiet_NaN();
    }

    return sum / static_cast<T>(count);
}

/*
 Funkcja wyszukująca indeks minimalnej wartości błędu.

 Minimalny błąd bezwzględny odpowiada granicy pomiędzy:
 - obszarem dominacji błędu obcięcia (dla większych h),
 - obszarem dominacji błędu zaokrągleń (dla bardzo małych h).

 W praktyce jest to punkt, od którego dalsze zmniejszanie
 kroku siatki nie poprawia dokładności obliczeń.
*/
template<typename T>
size_t find_min_error_index(const std::vector<T>& errors)
{
    size_t idx = 0;
    for (size_t i = 1; i < errors.size(); ++i)
        if (errors[i] < errors[idx])
            idx = i;
    return idx;
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
 Funkcja run_for_type realizuje pełny eksperyment numeryczny
 dla zadanego typu zmiennoprzecinkowego (double lub long double).

 Dla kolejnych kroków siatki h:
 - oblicza przybliżenia pochodnej w punktach:
   * x0 = 0        (początek przedziału),
   * xc = π/4      (środek przedziału),
   * xn = π/2      (koniec przedziału),
 - wyznacza błędy bezwzględne,
 - zapisuje wartości log10(h) oraz log10(|e|),
 - oblicza lokalne rzędy dokładności,
 - wyznacza uśrednione rzędy doświadczalne,
 - identyfikuje wpływ błędów maszynowych,
 - zapisuje dane do pliku do dalszej analizy graficznej.
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
    std::vector<T> p_f2, p_f3, p_c2, p_b2, p_b3;

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

        /*
            Obliczenie lokalnych rzędów dokładności następuje
            po zgromadzeniu co najmniej dwóch wartości błędu.

            Każdy rząd p_i odpowiada nachyleniu krzywej
            log|e(h)| względem log(h) pomiędzy dwoma kolejnymi krokami.
        */
        if(logh.size() >= 2)
        {
            size_t i = logh.size() - 1;
            T temp_accuracy_order = static_cast<T>(0);

            temp_accuracy_order = experimental_order(
                std::pow(static_cast<T>(10), ef2[i - 1]),
                std::pow(static_cast<T>(10), ef2[i]),
                std::pow(static_cast<T>(10), logh[i - 1]),
                std::pow(static_cast<T>(10), logh[i]));
            p_f2.push_back(temp_accuracy_order);

            temp_accuracy_order = experimental_order(
                std::pow(static_cast<T>(10), ef3[i - 1]),
                std::pow(static_cast<T>(10), ef3[i]),
                std::pow(static_cast<T>(10), logh[i - 1]),
                std::pow(static_cast<T>(10), logh[i]));
            p_f3.push_back(temp_accuracy_order);

            temp_accuracy_order = experimental_order(
                std::pow(static_cast<T>(10), ec2[i - 1]),
                std::pow(static_cast<T>(10), ec2[i]),
                std::pow(static_cast<T>(10), logh[i - 1]),
                std::pow(static_cast<T>(10), logh[i]));
            p_c2.push_back(temp_accuracy_order);

            temp_accuracy_order = experimental_order(
                std::pow(static_cast<T>(10), eb2[i - 1]),
                std::pow(static_cast<T>(10), eb2[i]),
                std::pow(static_cast<T>(10), logh[i - 1]),
                std::pow(static_cast<T>(10), logh[i]));
            p_b2.push_back(temp_accuracy_order);

            temp_accuracy_order = experimental_order(
                std::pow(static_cast<T>(10), eb3[i - 1]),
                std::pow(static_cast<T>(10), eb3[i]),
                std::pow(static_cast<T>(10), logh[i - 1]),
                std::pow(static_cast<T>(10), logh[i]));
            p_b3.push_back(temp_accuracy_order);
        }
    }

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n" << "=========================================\n"
            << "RZĘDY DOKŁADNOŚCI PRZYBLIŻEŃ RÓŹNICOWYCH\n"
            << "=========================================\n";

    std::cout << "=== Forward 2 ===\n"
            << "Rząd teoretyczny: O(h)\n"
            << "Rząd doświadczalny: "
            << average_experimental_order(p_f2, 5, 20) << "\n";
    std::cout << "=== Forward 3 ===\n"
            << "Rząd teoretyczny: O(h^2)\n"
            << "Rząd doświadczalny: "
            << average_experimental_order(p_f3, 5, 20) << "\n";
    std::cout << "=== Central 2 ===\n"
            << "Rząd teoretyczny: O(h^2)\n"
            << "Rząd doświadczalny: "
            << average_experimental_order(p_c2, 5, 20) << "\n";
    std::cout << "=== Backward 2 ===\n"
            << "Rząd teoretyczny: O(h)\n"
            << "Rząd doświadczalny: "
            << average_experimental_order(p_b2, 5, 20) << "\n";
    std::cout << "=== Backward 3 ===\n"
            << "Rząd teoretyczny: O(h^2)\n"
            << "Rząd doświadczalny: "
            << average_experimental_order(p_b3, 5, 20) << "\n";

    std::cout << std::scientific << std::setprecision(6);
    
    size_t idx_min_f2 = find_min_error_index(ef2);
    size_t idx_min_f3 = find_min_error_index(ef3);
    size_t idx_min_c2 = find_min_error_index(ec2);
    size_t idx_min_b2 = find_min_error_index(eb2);
    size_t idx_min_b3 = find_min_error_index(eb3);

    std::cout << "\n" << "=========================================\n"
            << "WPŁYW BŁĘDÓW MASZYNOWYCH\n"
            << "=========================================\n";

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Minimalny bład (forward 2) dla h ~ "
          << logh[idx_min_f2] << "\n";
    std::cout << "Minimalny bład (forward 3) dla h ~ "
          << logh[idx_min_f3] << "\n";
    std::cout << "Minimalny bład (central 2) dla h ~ "
          << logh[idx_min_c2] << "\n";
    std::cout << "Minimalny bład (backward 2) dla h ~ "
          << logh[idx_min_b2] << "\n";
    std::cout << "Minimalny bład (backward 3) dla h ~ "
          << logh[idx_min_b3] << "\n";
    
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