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

template<typename T>
int count_significant_digits(T value)
{
    // Convert number to scientific string with high precision
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(80) << value;
    std::string s = oss.str();

    // Separate mantissa from exponent
    auto pos_e = s.find('e');
    std::string mantissa = s.substr(0, pos_e);

    // Remove sign and decimal point → keep only digits
    std::string digits;
    for (char c : mantissa) {
        if (std::isdigit(c)) digits.push_back(c);
    }

    // 1. Remove leading zeros (numbers < 1)
    while (!digits.empty() && digits.front() == '0')
        digits.erase(digits.begin());

    // 2. Remove trailing zeros introduced by fixed precision
    while (!digits.empty() && digits.back() == '0')
        digits.pop_back();

    return digits.size();
}

template<typename T>
void analyze_type(const std::string &label)
{
    std::cout << "----- " << label << " -----\n";

    // experimental epsilon
    T eps_exp = compute_machine_epsilon_experimental<T>();

    // compute mantissa bits experimentally
    int mantissa_bits = -log2(eps_exp);

    // decimal digits estimated
    int decimal_digits_from_eps = count_significant_digits(eps_exp);

    std::cout << std::scientific << std::setprecision(70);
    std::cout << "experimental epsilon: " << static_cast<T>(eps_exp) << "\n";
    std::cout << "mantissa bits (experimental, including hidden bit if applicable): " << mantissa_bits << "\n";
    std::cout << "approx. decimal significant digits (from eps): " << decimal_digits_from_eps << "\n";
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