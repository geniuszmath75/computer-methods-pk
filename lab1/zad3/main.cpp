#include <iostream>
#include <cmath>
#include <iomanip>

const long double TOLX = 1e-15; // subsequent x approx. tolerance
const long double TOLF = 1e-15; // f(x) value tolerance
const int NMAX = 50; // max number of iteration

//*********************************************************
//* Functions for equation (a): fA(x) = tanh(x) + 2(x-1)
//*********************************************************
long double fA(long double x) {
    return std::tanh(x) + 2.0*(x - 1.0);
}

// Prime fA'(x) = (1 - tanh^2(x)) + 2
long double fAprime(long double x) {
    long double th = std::tanh(x);
    return (1.0 - th*th) + 2.0;
}

// Function gA(x) = 1 - 0.5*tanh(x), for Picard iteration
long double gA(long double x) {
    return 1.0 - 0.5*std::tanh(x);
}

//*********************************************************
//* Functions for equation (b): fA(x) = sinh(x) + x/4 - 1
//*********************************************************
long double fB(long double x) {
    return std::sinh(x) + (x / 4.0) - 1.0;
}

// Prime fB'(x) = cosh(x) + 1/4
long double fBprime(long double x) {
    return std::cosh(x) + 0.25;
}

// Function gB(x) = asinh(1 - x/4), for Picard iteration
long double gB(long double x) {
    return asinh(1.0 - x / 4.0);
}

//***************************************
//* 1. Picard iteration: x_{n+1} = g(x_n)
//***************************************
long double picard_iteration(
    long double (*f)(long double),
    long double (*g)(long double),
    long double (x0),
    const std::string &desc
)
{
    std::cout << "=== Picard iteration for " << desc << " ===\n";
    long double x = x0;

    for(int i = 1; i <= NMAX; i++)
    {
        long double fx = f(x);
        long double x_next = g(x);

        // Current iteration information
        std::cout << "Iteration " << i
        << ": x = " << std::setw(12) << x
        << ", f(x) = " << std::setw(12) << fx << "\n";

        // Checking stop conditions
        if(std::fabs(fx) < TOLF)
        {
            std::cout << "Specified precision has been reached f(x) < TOLF.\n\n";
            return x;
        }
        if(std::fabs(x_next - x) < TOLX)
        {
            std::cout << "Specified precision has been reached |x_{n+1} - x_n}| < TOLX.\n\n";
            return x_next;
        }

        x = x_next;
    }
    std::cout << "Maximum number of iteration reached in Picard iteration.\n\n";
    return x;
}

//***************************************
//* 2. Bisection method
//*    Requires [a, b] range, f(a)*f(b) < 0
//***************************************
long double bisection_method(
    long double (*f)(long double),
    long double a,
    long double b,
    const std::string &desc
)
{
    std::cout << "=== Bisection method for " << desc << " ===\n";

    long double fa = f(a);
    long double fb = f(b);

    if(fa * fb >= 0.0L)
    {
        std::cerr << "ERROR: f(a) i f(b) must have opposite signs (f(a)=" << fa << ", f(b)=" << fb << ").\n";
        return 0.0L;
    }

    for(int i = 1; i <= NMAX; i++)
    {
        long double m = 0.5 *(a + b);
        long double fm = f(m);

        // Current iteration information
        std::cout << "Iteration " << i 
        << ": a = " << std::setw(12) << a 
        << ", b = " << std::setw(12) << b 
        << ", m = " << std::setw(12) << m 
        << ", f(m) = " << std::setw(12) << fm << "\n";

        // Checking stop conditions
        if(std::fabs(fm) < TOLF)
        {
            std::cout << "Specified precision has been reached f(m) < TOLF.\n\n";
            return m;
        }
        if(std::fabs(b - a) < TOLX)
        {
            std::cout << "Specified precision has been reached (b - a) < TOLX.\n\n";
            return m;
        }

        if(fa * fm < 0.0)
        {
            b = m;
            fb = fm;
        }
        else
        {
            a = m;
            fa = fm;
        }
    }
    std::cout << "Maximum number of iteration reached in Bisection method.\n\n";
    return 0.5*(a + b);
}

//**************************************************
//* 3. Newton method: x_{n+1} = x_n - f(x_n)/f'(x_n)
//**************************************************
long double newton_method(
    long double (*f)(long double),
    long double (*fprime)(long double),
    long double x0,
    const std::string &desc
) {
    std::cout << "=== Newton method for " << desc << " ===\n";
    long double x = x0;

    for(int i = 1; i <= NMAX; i++) {
        long double fx = f(x);
        long double fpx = fprime(x);

        std::cout << "Iteration " << i 
                  << ": x = " << std::setw(12) << x 
                  << ", f(x) = " << std::setw(12) << fx 
                  << ", f'(x) = " << std::setw(12) << fpx << "\n";

        if (std::fabs(fpx) < 1e-14) {
            std::cerr << "Derivative close to zero - calculation cancelled.\n\n";
            return x;
        }
        if (std::fabs(fx) < TOLF) {
            std::cout << "Specified precision has been reached f(x) < TOLF.\n\n";
            return x;
        }

        long double x_next = x - fx / fpx;

        if (std::fabs(x_next - x) < TOLX) {
            std::cout << "Specified precision has been reached |x_{n+1} - x_n| < TOLX.\n\n";
            return x_next;
        }
        x = x_next;
    }
    std::cout << "Maximum number of iteration reached in Newton method.\n\n";
    return x;
}

//*******************************************************************
//* 4. Secant method:
//*    x_{n+1} = x_n - f(x_n)*(x_n - x_{n-1}) / [f(x_n) - f(x_{n-1})]
//*******************************************************************
long double secant_method(
    long double (*f)(long double),
    long double x0,
    long double x1,
    const std::string &desc
) {
    std::cout << "=== Secant method for " << desc << " ===\n";
    
    long double f0 = f(x0);
    long double f1 = f(x1);

    for(int i = 1; i <= NMAX; i++) {
        if(std::fabs(f1 - f0) < 1e-14)
        {
            std::cerr << "f(x1) - f(x0) close to zero - calculation cancelled.";
            return x1;
        }

        std::cout << "Iteration " << i 
                  << ": x0 = " << std::setw(12) << x0
                  << ": x1 = " << std::setw(12) << x1 
                  << ", f(x0) = " << std::setw(12) << f0 
                  << ", f(x1) = " << std::setw(12) << f1 << "\n";

        long double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        long double f2 = f(x2);

        if (std::fabs(f2) < TOLF) {
            std::cout << "Specified precision has been reached f(x2) < TOLF.\n\n";
            return x2;
        }
        if (std::fabs(x2 - x1) < TOLX) {
            std::cout << "Specified precision has been reached |x_{n+1} - x_n| < TOLX.\n\n";
            return x2;
        }
        
        x0 = x1;
        x1 = x2;
        f0 = f1;
        f1 = f2;
    }

    std::cout << "Maximum number of iteration reached in Secant method.\n\n";
    return x1;
}

int main()
{
    std::cout << std::fixed << std::setprecision(15);

    //***************************************
    //* Equation: (a): tanh(x) + 2(x - 1) = 0
    //***************************************
    std::cout << "\n******** Equation: (a): tanh(x) + 2(x - 1) = 0 ********\n";

    // 1. Picard: x_{n+1} = 1 - 0.5*tanh(x_n)
    {
        long double x0 = 0.5L;
        long double root = picard_iteration(fA, gA, x0, "(a)");
        std::cout << "Picard solution (a) ~ " << root << "\n\n";
    }
    
    // 2. Bisection on the interval [0, 1]: 
    {
        long double a = 0.0L, b = 1.0L;
        long double root = bisection_method(fA, a, b, "(a)");
        std::cout << "Bisection solution (a) ~ " << root << "\n\n";
    }

    // 3. Newton - start x0 = 0.5
    {
        long double x0 = 0.5L;
        long double root = newton_method(fA, fAprime, x0, "(a)");
        std::cout << "Newton solution (a) ~ " << root << "\n\n";
    }

    // 4. Secant - start x0 = 0.0, x1 = 1.0
    {
        long double x0 = 0.0L;
        long double x1 = 1.0L;
        long double root = secant_method(fA, x0, x1, "(a)");
        std::cout << "Secant solution (a) ~ " << root << "\n\n";
    }

    //***************************************
    //* Equation: (b): sinh(x) + (x/4) - 1 = 0
    //***************************************
    std::cout << "\n******** Equation: (b): sinh(x) + (x/4) - 1 = 0 ********\n";

    // 1. Picard: x_{n+1} = asinh(1 - x_n/4)
    {
        long double x0 = 1.0L;
        long double root = picard_iteration(fB, gB, x0, "(b)");
        std::cout << "Picard solution (b) ~ " << root << "\n\n";
    }
    
    // 2. Bisection on the interval [0, 2]:
    //    fB(0)=-1, fB(2)=3.1269 => sign change 
    {
        long double a = 0.0L, b = 2.0L;
        long double root = bisection_method(fB, a, b, "(b)");
        std::cout << "Bisection solution (b) ~ " << root << "\n\n";
    }

    // 3. Newton - start x0 = 1.0
    {
        long double x0 = 1.0L;
        long double root = newton_method(fB, fBprime, x0, "(b)");
        std::cout << "Newton solution (b) ~ " << root << "\n\n";
    }

    // 4. Secant - start x0 = 0.0, x1 = 2.0
    {
        long double x0 = 0.0L;
        long double x1 = 2.0L;
        long double root = secant_method(fB, x0, x1, "(b)");
        std::cout << "Secant solution (b) ~ " << root << "\n\n";
    }

    return 0;
}