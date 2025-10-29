#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#define USE_STABLE_METHOD

static const long double EPSILON = 1.084e-19L;

struct DataPoint
{
    long double log10x;
    long double x;
    long double f_exact;
};


template <typename T>
T compute_unstable_aproximate_function_result(T x)
{
    T six = static_cast<T>(6);
    return x * x * x / (six * (std::sinh(x) - x));
}

template <typename T>
T  compute_stable_aproximate_function_result(T x)
{
    T six = static_cast<T>(6);
    T zero = static_cast<T>(0);
    T two = static_cast<T>(2);
    T three = static_cast<T>(3);
    
    T term = x * x * x / six;
    T sum = zero;

    for(int n = 1; n <= 10; ++n)
    {
        sum += term;
        term *= (x * x) / ((two * n + two) * (two * n + three));
    }
    return x * x * x / (six * sum);
}

int main()
{
    const std::string inputFile = "dane_do_laboratorium_2.txt";
    const std::string outputFile = "relative_errors.txt";

    std::ifstream fin(inputFile);
    if(!fin.is_open())
    {
        std::cerr << "Error: cannot open input file: " << inputFile << "\n";
        return 1;
    }

    std::ofstream fout(outputFile);
    if(!fout.is_open())
    {
        std::cerr << "Error: cannot open output file: " << outputFile << "\n";
        return 1;
    }

    fout << "====================================================================\n";
    fout << "log10(x)\tlog10(|relative_error|)\n";
    fout << "====================================================================\n";

    std::string line;
    while(std::getline(fin, line))
    {
        if(line.find('=') != std::string::npos || line.empty()) continue;

        DataPoint point;
        if(!(std::istringstream(line) >> point.log10x >> point.x >> point.f_exact)) continue;

        long double f_approx;

        #ifdef USE_STABLE_METHOD
            if(point.x < 0.9L)
            {
                f_approx = compute_stable_aproximate_function_result(point.x);    
            }
            else
            {
                f_approx = compute_unstable_aproximate_function_result(point.x);
            }
        #else
            f_approx = compute_unstable_aproximate_function_result(point.x);
        #endif
        

        // Calculate relative error
        long double relativeError;
        if(f_approx == point.f_exact)
        {
            relativeError = EPSILON;
        } else
        {
            relativeError = (fabsl(f_approx - point.f_exact) / fabsl(point.f_exact));
        }

        fout << std::fixed << std::setprecision(5) << point.log10x
             << "\t" << std::setprecision(8)
             << std::log10(relativeError) << "\n";
    }

    fin.close();
    fout.close();

    std::cout << "Results saved to " << outputFile << "\n";

    return 0;
}