#include "ode-suite.hpp"

ODE_Func_Long ode_fcn_long(User_Math_Long custom_equation)
{
    return [eq = std::move(custom_equation)](long double t, const std::vector<long double>& y)
    {
        if (y.empty())
        {
            throw std::invalid_argument("ode_fcn_long: State vector is empty");
        }

        size_t n = y.size();
        std::vector<long double> dydt(n);

        for (size_t i = 0; i + 1 < n; ++i)
        {
            dydt[i] = y[i + 1];
        }

        dydt[n - 1] = eq(t, y);

        return dydt;
    };
}