#include "ode-suite.hpp"


ODE_Func ode_fcn(User_Math custom_equation)
{
	return [eq = std::move(custom_equation)](double t, const std::vector<double>& y)
	{
		if (y.empty())
		{
			throw std::invalid_argument("ode_fcn: State vector is empty");
		}

		size_t n = y.size();
		std::vector<double> dydt(n);

		for (size_t i = 0; i + 1 < n; ++i)
		{
			dydt[i] = y[i + 1];
		}

		dydt[n - 1] = eq(t, y);

		return dydt;
	};
}