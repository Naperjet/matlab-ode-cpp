#include "ode-suite.hpp"
#include <string>

namespace
{
	std::vector<double> bs23_hermite(double s, double h, const std::vector<double>& y, const std::vector<double>& y_next, const std::vector<double>& k0, const std::vector<double>& k3)
	{
		double s2 = s * s;
		double h00 = (1.0 - s) * (1.0 - s) * (1.0 + 2.0 * s);
		double h10 = s * (1.0 - s) * (1.0 - s);
		double h01 = s2 * (3.0 - 2.0 * s);
		double h11 = s2 * (s - 1.0);

		size_t n = y.size();
		std::vector<double> interp(n);
		for (size_t i = 0; i < n; ++i)
		{
			interp[i] = h00 * y[i] + h01 * y_next[i] + h * (h10 * k0[i] + h11 * k3[i]);
		}

		return interp;
	}

	double estimate_initial_step(double t0, double tf, const std::vector<double>& y0, const std::vector<double>& f0, double reltol, double abstol)
	{
		size_t n = y0.size();
		double norm_y = 0.0;
		double norm_f = 0.0;

		for (size_t i = 0; i < n; ++i)
		{
			double sc = std::max(abstol, reltol * std::abs(y0[i]));
			norm_y = std::max(norm_y, std::abs(y0[i]) / sc);
			norm_f = std::max(norm_f, std::abs(f0[i]) / sc);
		}

		double h0 = (norm_y < 1e-5 || norm_f < 1e-5) ? 1e-6 : 0.01 * (norm_y / norm_f);

		return std::min(h0, std::abs(tf - t0));
	}
}

ODE_Res ode23(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options = ODE_Set{})
{
	if (tspan.size() < 2)
	{
		throw std::invalid_argument("ode23: tspan must contain at least two elements");
	}

	bool increasing = (tspan[1] >= tspan[0]);
	for (size_t i = 1; i < tspan.size(); ++i)
	{
		if (increasing && tspan[i] < tspan[i - 1])
		{
			throw std::invalid_argument("ode23: tspan is not strictly increasing");
		}

		if (!increasing && tspan[i] > tspan[i - 1])
		{
			throw std::invalid_argument("ode23: tspan is not strictly decreasing");
		}
	}

	ODE_Res res;
	const bool return_specific_points = (tspan.size() > 2);
	const double t0 = tspan.front();
	const double tf = tspan.back();

	double t = t0;
	std::vector<double> y = Init_Value;
	const size_t n = y.size();

	if (n == 0)
	{
		throw std::invalid_argument("ode23: Initial state cannot be empty.");
	}

	res.Time.push_back(t);
	res.Result.push_back(y);

	double direction = (tf >= t0) ? 1.0 : -1.0;

	double h;
	if (Options.Init_Step > 0.0)
	{
		h = Options.Init_Step;
	}
	else
	{
		auto f0 = F(t, y);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}
		h = estimate_initial_step(t0, tf, y, f0, Options.Rel_Tol, Options.Abs_Tol);
	}

	h = direction * std::abs(h);

	double h_max = (Options.Max_Step > 0.0) ? Options.Max_Step : std::abs(tf - t0) / 10.0;
	h_max = std::abs(h_max);

	const double eps = std::numeric_limits<double>::epsilon();
	const size_t max_iter = 1'000'000;

	std::vector<std::vector<double>> k(4, std::vector<double>(n));
	std::vector<double> y_next(n);

	k[0] = F(t, y);
	if (Options.Stats)
	{
		res.stats.nfevals++;
	}

	size_t tspan_idx = 1;
	size_t iter = 0;

	while ((direction > 0.0 && t < tf) || (direction<0.0 && t>tf))
	{
		if (++iter > max_iter)
		{
			throw std::runtime_error("ode23: Consider stiff solvers. Maximum iteration count exceeded.");
		}

		if ((direction > 0.0 && t + h > tf) || (direction < 0.0 && t + h < tf))
		{
			h = tf - t;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * (1.0 / 2.0) * k[0][i];
		}
		k[1] = F(t + h * (1.0 / 2.0), y_next);
		if (Options.Stats) 
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * (3.0 / 4.0) * k[1][i];
		}
		k[2] = F(t + h * (3.0 / 4.0), y_next);
		if (Options.Stats) 
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * ((2.0 / 9.0) * k[0][i] + (1.0 / 3.0) * k[1][i] + (4.0 / 9.0) * k[2][i]);
		}
		k[3] = F(t + h, y_next);
		if (Options.Stats) 
		{
			res.stats.nfevals++;
		}

		double max_err_ratio = 0.0;
		for (size_t i = 0; i < n; ++i)
		{
			double error = h * ((-5.0 / 72.0) * k[0][i] + (1.0 / 12.0) * k[1][i] + (1.0 / 9.0) * k[2][i] - (1.0 / 8.0) * k[3][i]);
			double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y[i]), std::abs(y_next[i])));

			max_err_ratio = std::max(max_err_ratio, std::abs(error) / tol);

			if (!std::isfinite(y_next[i]))
			{
				throw std::runtime_error(
					"ode23: solution became non-finite at t = " + std::to_string(t + h));
			}
		}

		if (max_err_ratio <= 1.0)
		{
			if (Options.Stats)
			{
				res.stats.nsteps++;
			}

			if (!Options.NonNegative.empty())
			{
				for (int idx : Options.NonNegative)
				{
					if (idx >= 0 && idx < static_cast<int>(n))
					{
						y_next[idx] = std::max(y_next[idx], 0.0);
					}
				}
			}

			if (return_specific_points)
			{
				while (tspan_idx < tspan.size() && ((direction > 0.0 && tspan[tspan_idx] <= t + h) || (direction < 0.0 && tspan[tspan_idx] >= t + h)))
				{
					double s = (tspan[tspan_idx] - t) / h;
					res.Time.push_back(tspan[tspan_idx]);
					res.Result.push_back(bs23_hermite(s, h, y, y_next, k[0], k[3]));
					++tspan_idx;
				}
			}

			else
			{
				if (Options.Refine > 1)
				{
					for (int j = 1; j < Options.Refine; ++j)
					{
						double s = static_cast<double>(j) / static_cast<double>(Options.Refine);
						res.Time.push_back(t + s * h);
						res.Result.push_back(bs23_hermite(s, h, y, y_next, k[0], k[3]));
					}
				}

				res.Time.push_back(t + h);
				res.Result.push_back(y_next);
			}

			if (Options.Events)
			{
				auto ev_prev = (*Options.Events)(t, y);
				auto ev_next = (*Options.Events)(t + h, y_next);

				for (size_t ev_idx = 0; ev_idx < ev_prev.value.size(); ++ev_idx)
				{
					double v0 = ev_prev.value[ev_idx];
					double v1 = ev_next.value[ev_idx];

					bool crossed = false;
					if (v0 * v1 < 0.0)
					{
						int dir = (ev_idx < ev_prev.direction.size()) ? ev_prev.direction[ev_idx] : 0;
						if (dir == 0)
						{
							crossed = true;
						}
						else if (dir > 0 && v1 > v0) 
						{
							crossed = true;
						}
						else if (dir < 0 && v1 < v0) 
						{
							crossed = true;
						}
					}

					if (crossed)
					{
						double s_lo = 0.0, s_hi = 1.0;
						double v_lo = v0, v_hi = v1;

						for (int b = 0; b < 30; ++b)
						{
							double s_mid = 0.5 * (s_lo + s_hi);
							auto y_mid = bs23_hermite(s_mid, h, y, y_next, k[0], k[3]);
							auto ev_mid = (*Options.Events)(t + s_mid * h, y_mid);
							double v_mid = ev_mid.value[ev_idx];

							if (v_lo * v_mid <= 0.0)
							{
								s_hi = s_mid;
								v_hi = v_mid;
							}
							else
							{
								s_lo = s_mid;
								v_lo = v_mid;
							}
						}

						double s_ev = 0.5 * (s_lo + s_hi);
						auto y_ev = bs23_hermite(s_ev, h, y, y_next, k[0], k[3]);
						double t_ev = t + s_ev * h;

						res.Te.push_back(t_ev);
						res.Ye.push_back(y_ev);
						res.Ie.push_back(static_cast<int>(ev_idx) + 1);

						bool terminal = (ev_idx < ev_prev.isterminal.size()) ? ev_prev.isterminal[ev_idx] : false;
						if (terminal)
						{
							t = t_ev;
							y = y_ev;
							return res;
						}
					}
				}
			}

			t += h;
			y = y_next;
			k[0] = k[3];

			if (return_specific_points && tspan_idx >= tspan.size())
			{
				break;
			}
		}

		else
		{
			if (Options.Stats)
			{
				res.stats.nfailed++;
			}
		}

		double h_scale = 0.9 * std::pow(1.0 / std::max(max_err_ratio, 1e-10), 1.0 / 3.0);
		h_scale = std::clamp(h_scale, 0.1, 5.0);

		h = direction * std::max(std::abs(h * h_scale), 16.0 * eps * std::abs(t));
		h = direction * std::min(std::abs(h), h_max);

		if (std::abs(h) < 16.0 * eps * std::abs(t) && max_err_ratio > 1.0)
		{
			throw std::runtime_error("ode23: Step Size below minimum. Consider a stiff solver");
		}
	}

	return res;
}