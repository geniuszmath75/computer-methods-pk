#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#define USE_STABLE_METHOD

struct DataPoint
{
    long double log10x;
    long double x;
    long double f_exact;
};


template <typename T>
T compute_unstable_aproximate_function_result(T x)
{
    T numerator = x * x * x;
    T denominator = 6.0 * (std::sinh(x) - x);
    return numerator / denominator;
}

template <typename T>
T  compute_stable_aproximate_function_result(T x)
{
    if(std::abs(x) < 1e-3)
    {
        T x2 = x * x;
        T x4 = x2 * x2;
        return 1.0 - (x2 / 20.0) + (x4 / 840.0);
    }
    else
    {
        T numerator = x * x * x;
        T denominator = 6.0 * (std::sinh(x) - x);
        return numerator / denominator;
    }
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
    fout << "log10(x)\tlog(|relative_error|)\trelative_error\tf_exact\tf_approx\n";
    fout << "====================================================================\n";

    std::string line;
    while(std::getline(fin, line))
    {
        if(line.find('=') != std::string::npos || line.empty()) continue;

        DataPoint point;
        if(!(std::istringstream(line) >> point.log10x >> point.x >> point.f_exact)) continue;

        #ifdef USE_STABLE_METHOD
            long double f_approx = compute_stable_aproximate_function_result(point.x);    
        #else
            long double f_approx = compute_unstable_aproximate_function_result(point.x);
        #endif
        

        // Calculate relative error
        long double relativeError = std::abs((f_approx - point.f_exact) / point.f_exact);

        fout << std::fixed << std::setprecision(5) << point.log10x
             << "\t" << std::scientific << std::setprecision(20)
             << std::log10(relativeError) << "\t"
             << relativeError << "\t"
             << point.f_exact << "\t"
             << f_approx << "\n"; 
    }

    fin.close();
    fout.close();

    std::cout << "Results saved to " << outputFile << "\n";

    return 0;
}