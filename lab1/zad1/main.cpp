#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <type_traits>
#include <string>

// Compute machine epsilon experimentally for type T
template<typename T>
T compute_machine_epsilon_experimental() 
{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type");
    // Start from 1 and repeatedly halve until 1 + eps/2 == 1
    T eps = static_cast<T>(1);
    T one = static_cast<T>(1);
    while ((one + eps / static_cast<T>(2)) != one)
    {
        eps = eps / static_cast<T>(2);
    }
    return eps;
    
}

// Compute mantissa bits from epsilon (eps = 2^(1 - p) => p = 1 - log2(eps))
template <typename T>
T compute_mantissa_bits_from_epsilon(T eps) {
    // compute log2 from epsilon
    long double log2_eps = std::log2(static_cast<long double>(eps));
    long double p = static_cast<long double>(1) - log2_eps;

    // round to nearest int
    int p_int = static_cast<int>(std::llround(p));
    return p_int;
}

// Compute appproximate decimal significant digits from epsilon
template <typename T>
T compute_decimal_significant_digits_from_epsilon(T eps) {
    long double neg_log10 = -std::log10(static_cast<long double>(eps));

    // round down to nearest whole digit as a conservative estimate
    int digits = static_cast<int>(std::floor(neg_log10 + 1e-12L));
    return digits;
}

template<typename T>
void analyze_type(const std::string &label)
{
    std::cout << "----- " << label << " -----\n";

    // experimental epsilon
    T eps_exp = compute_machine_epsilon_experimental<T>();

    // library epsilon for comparison
    long double eps_std = static_cast<long double>(std::numeric_limits<T>::epsilon());

    // compute mantissa bits experimentally
    int mantissa_bits = compute_mantissa_bits_from_epsilon(eps_exp);

    // decimal digits estimated
    int decimal_digits_from_eps = compute_decimal_significant_digits_from_epsilon(eps_exp);

    // alternative estimate using mantissa bits
    long double digits_est_from_bits = mantissa_bits * std::log10(static_cast<long double>(2));

    std::cout << std::scientific << std::setprecision(10);
    std::cout << "experimental epsilon: " << static_cast<long double>(eps_exp) << "\n";
    std::cout << "std::numeric_limits<" << label << ">::epsilon(): " << eps_std << "\n";
    std::cout << std::defaultfloat;
    std::cout << "mantissa bits (experimental, including hidden bit if applicable): " << mantissa_bits << "\n";
    std::cout << "approx. decimal significant digits (from eps): " << decimal_digits_from_eps << "\n";
    std::cout << "approx. decimal digits (from mantissa bits): " << std::fixed << std::setprecision(3) << digits_est_from_bits << "\n";

    // Show relationship check: verify that 1 + eps == nextafter(1, 2)
    T one = static_cast<T>(1);
    T next_after = std::nextafter(one, static_cast<T>(2)); // next representable > 1
    T diff = next_after - one;
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "nextafter(1, 2) - 1 = " << static_cast<long double>(diff) << "\n";
    std::cout << "Is experimental epsilon equal to nextafter(1, 2) - 1? " << ( std::abs(static_cast<long double>(eps_exp - diff)) <= std::numeric_limits<long double>::epsilon() * std::fabs(static_cast<long double>(diff) + static_cast<long double>(eps_exp)) ? "YES" : "NO" ) << "\n";

    std::cout << std::defaultfloat << std::setprecision(6) << "\n";
}

int main()
{
    std::cout << "Machine epsilon experimental discovery \n\n";

    // Analyze float
    analyze_type<float>("float");

    // Analyze double
    analyze_type<double>("double");

    // Analyze long double
    analyze_type<long double>("long double");

    return 0;
}