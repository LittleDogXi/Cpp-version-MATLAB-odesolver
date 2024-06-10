#include"odesolver.hpp"
#include<stdio.h>
using std::cout;
using std::endl;
using namespace odesolver;

void odefun2(double t, const state_type& x, state_type& dxdt, void* auxdata);

typedef long double ld;

int main()
{
    state_type y0({ 0.1, 0.2 });
    //state_type sol;
    OdeSol sol;

    double atol = 1.0e-9, rtol = atol;
    OdeSet options(atol, rtol);

    ode45<2>(odefun2, 0.0, 5.0, y0, &sol, &options, nullptr);
    cout << sol.tout.size() << endl;
    printf("%.20e\n", sol.yout.back()[0]);
    printf("%.20e\n", sol.yout.back()[1]);

    return 0;
}

void odefun2(double t, const state_type& x, state_type& dxdt, void* auxdata)
{
    double x1 = x[0], x2 = x[1];
    dxdt[0] = x1 + x2 * sin(t);
    dxdt[1] = x1 - x2 * exp(-t);
}
