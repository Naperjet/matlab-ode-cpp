#include "ode-suite.hpp"

ODEFunc ode_fcn(UserMath custom_equation)
{
    ODEFunc F = [custom_equation](double t, const std::vector<double>& y)
    {
        size_t n = y.size();

        std::vector<double> dydt(n);

        for (int i = 0; i < n - 1; i++)
        {
            dydt[i] = y[i + 1];
        }

        dydt[n - 1] = custom_equation(t,y);

        return dydt;
    };

    return F;
}