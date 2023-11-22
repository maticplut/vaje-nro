#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

double exp3x(double x) {
    return std::exp(3 * x);
}

double arctanTaylor(double x, int terms) {
    double rezultat = 0.0;
    for (int n = 0; n < terms; ++n) {
        rezultat += std::pow(-1, n) / (2 * n + 1) * std::pow(x / 2, 2 * n + 1);
    }
    return rezultat;
}

double trapezoidalIntegration(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double vsota = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        vsota += f(a + i * h);
    }
    return h * vsota;
}

int main() {
    
    double spodnja = 0.0;
    double zgornja = M_PI / 4.0;

    int n;
    std::cout << "Stevilo podintervalov (n): "; 
    std::cin >> n;

    double rezultat = trapezoidalIntegration([](double x) {
        return exp3x(x) * arctanTaylor(x, 20); //20 členov
        }, spodnja, zgornja, n);

    std::cout << "Rezultat trapezne integracije za funkcijo e^(3x)*arctan(x/2): " << rezultat << std::endl;

    return 0;
}