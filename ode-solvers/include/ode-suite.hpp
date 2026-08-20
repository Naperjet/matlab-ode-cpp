#pragma once

#include <vector>
#include <functional>
#include <optional>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>

using ODE_Func = std::function<std::vector<double>(double, const std::vector<double>&)>;

using User_Math = std::function<double(double, const std::vector<double>&)>;

struct Event_Result
{
	std::vector<double> value;
	std::vector<bool> isterminal;
	std::vector<int> direction;
};

using Event_Func = std::function<Event_Result(double, const std::vector<double>&)>;

struct ODE_Res
{
	std::vector<double> Time;
	std::vector<std::vector<double>> Result;

	std::vector<double> Te;
	std::vector<std::vector<double>> Ye;
	std::vector<int> Ie;

	struct Stats
	{
		int nsteps = 0;
		int nfailed = 0;
		int nfevals = 0;
	} stats;
};

struct ODE_Set
{
	double Rel_Tol = 1e-3;
	double Abs_Tol = 1e-6;

	double Init_Step = 0.0;
	double Max_Step = 0.0;

	int Refine = 4;

	std::optional<Event_Func> Events;

	std::vector<int> NonNegative;

	bool Stats = false;
};

ODE_Func ode_fcn(User_Math custom_equation);

ODE_Res ode45(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode23(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode15s(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode23s(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode23t(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode23tb(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode113(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});
ODE_Res ode78(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{});

using ODE_Func_Long = std::function<std::vector<long double>(long double, const std::vector<long double>&)>;

using User_Math_Long = std::function<long double(long double, const std::vector<long double>&)>;

struct Event_Result_Long
{
	std::vector<long double> value;
	std::vector<bool> isterminal;
	std::vector<int> direction;
};

using Event_Func_Long = std::function<Event_Result_Long(long double, const std::vector<long double>&)>;

struct ODE_Res_Long
{
	std::vector<long double> Time;
	std::vector<std::vector<long double>> Result;

	std::vector<long double> Te;
	std::vector<std::vector<long double>> Ye;
	std::vector<int> Ie;

	struct Stats
	{
		int nsteps = 0;
		int nfailed = 0;
		int nfevals = 0;
	} stats;
};

struct ODE_Set_Long
{
	long double Rel_Tol = 1e-9L;
	long double Abs_Tol = 1e-16L;

	long double Init_Step = 0.0L;
	long double Max_Step = 0.0L;

	int Refine = 4;

	std::optional<Event_Func_Long> Events;

	std::vector<int> NonNegative;

	bool Stats = false;
};

ODE_Func_Long ode_fcn_long(User_Math_Long custom_equation);
ODE_Res_Long ode113_long(const ODE_Func_Long& F, const std::vector<long double>& tspan, const std::vector<long double>& Init_Value, const ODE_Set_Long& Options = ODE_Set_Long{});
ODE_Res ode89(const ODE_Func_Long& F, const std::vector<long double>& tspan, const std::vector<long double>& Init_Value, const ODE_Set_Long& Options = ODE_Set_Long{});
