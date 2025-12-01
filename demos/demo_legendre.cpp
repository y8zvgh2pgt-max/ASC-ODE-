#include <iostream>
#include <fstream>
#include <cmath>

#include <autodiff.hpp>

using namespace std;
using namespace ASC_ode;

template <typename T>
T Legendre(int n, T x) {
    if (n == 0) return T(1.0);
    if (n == 1) return x;
    
    T p_prev2 = 1.0;
    T p_prev1 = x;
    T p_curr = x;

    for (int i = 2; i <= n; ++i) {
        p_curr = (T(2 * i - 1) * x * p_prev1 - T(i - 1) * p_prev2) / T(i);
        p_prev2 = p_prev1;
        p_prev1 = p_curr;
    }
    return p_curr;
}

int main() {
    ofstream out("output_legendre.txt");
    
    int max_order = 5;
    int steps = 200;
    double x_start = -1.0;
    double x_end = 1.0;
    double dx = (x_end - x_start) / (steps - 1);

    out << "x";
    for (int n = 0; n <= max_order; n++) out << " P" << n << " dP" << n;
    out << endl;

    for (int i = 0; i < steps; i++) {
        double x_val = x_start + i * dx;
        
        AutoDiff<1> x = Variable<0>(x_val);

        out << x_val;

        for (int n = 0; n <= max_order; n++) {
            AutoDiff<1> pn = Legendre(n, x);
            
            out << " " << pn.value() << " " << pn.deriv()[0];
        }
        out << endl;
    }

    return 0;
}