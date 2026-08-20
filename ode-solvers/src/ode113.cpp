#include "ode-suite.hpp"
#include <string>

namespace
{
	constexpr double G_Coeff[9] =
	{
		1.0,
		1.0 / 2.0,
		5.0 / 12.0,
		3.0 / 8.0,
		251.0 / 720.0,
		95.0 / 288.0,
		19087.0 / 60480.0,
		5257.0 / 17280.0,
		1070017.0 / 3628800.0
	};

	std::vector<double> hermite_interp(double s, double h, const std::vector<double>& y, const std::vector<double>& f, const std::vector<double>& y_next, const std::vector<double>& f_next)
	{
		const double s2 = s * s;
		const double s3 = s2 * s;
		const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
		const double h10 = s3 - 2.0 * s2 + s;
		const double h01 = -2.0 * s3 + 3.0 * s2;
		const double h11 = s3 - s2;

		size_t n = y.size();
		std::vector<double> interp(n);
		for (size_t i = 0; i < n; ++i)
		{
			interp[i] = h00 * y[i] + h01 * y_next[i] + h * (h10 * f[i] + h11 * f_next[i]);
		}

		return interp;
	}

	double estimate_initial_step(double t0, double tf, const std::vector<double>& y0, const std::vector<double>& f0, double Reltol, double AbsTol)
	{
		size_t n = y0.size();
		double norm_y = 0.0;
		double norm_f = 0.0;

		for (int i = 0; i < n; ++i)
		{
			double sc = std::max(AbsTol, Reltol * std::abs(y0[i]));
			norm_y = std::max(norm_y, std::abs(y0[i]) / sc);
			norm_f = std::max(norm_f, std::abs(f0[i]) / sc);
		}

		double h0 = (norm_y < 1e-5 || norm_f < 1e-5) ? 1e-6 : 0.01 * (norm_y / norm_f);
		return std::min(h0, std::abs(tf - t0));
	}
}

ODE_Res ode113(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options)
{
	if (tspan.size() < 2)
	{
		throw std::invalid_argument("ode113: tspan must contain at least two points");
	}

	bool increasing = (tspan[1] > tspan[0]);
	for (size_t i = 0; i < tspan.size(); ++i)
	{
		if (increasing && tspan[i] < tspan[i - 1])
		{
			throw std::invalid_argument("ode113: tspan must be strictly increasing.");
		}

		if (!increasing && tspan[i] > tspan[i - 1])
		{
			throw std::invalid_argument("ode113: tspan must be strictly decreasing.");
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
		throw std::invalid_argument("ode113: initial state cannot be empty");
	}

	res.Time.push_back(t);
	res.Result.push_back(y);

	double dir = (tf > t0) ? 1.0 : -1.0;

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

	h = dir * std::abs(h);

	double h_max = (Options.Max_Step > 0.0) ? Options.Max_Step : std::abs(tf - t0) / 10;
	h_max = std::abs(h_max);

	const double eps = std::numeric_limits<double>::epsilon();
	const size_t max_iter = 1'000'000;

	const int max_order = 7;

	std::vector<std::vector<double>> phi(max_order + 2, std::vector<double>(n, 0.0));

	std::vector<double> f_now = F(t, y);
	if (Options.Stats)
	{
		res.stats.nfevals++;
	}
	phi[0] = f_now;

	int order = 1;
	double h_old = h;
	size_t tspan_idx = 1;
	size_t iter = 0;
	bool done = false;

	while (!done)
	{
		if (++iter > max_iter)
		{
			throw std::runtime_error("ode113: maximum iterations exceeded.");
		}

		if ((dir > 0.0 && t + h > tf) || (dir < 0.0 && t + h < tf))
		{
			h = tf - t;
			done = true;
		}

		bool restarted = false;
		double h_ratio = std::abs(h / h_old);

		if ((h_ratio < 0.5 || h_ratio > 2.0) && order > 1)
		{
			order = 1;
			restarted = true;
			phi[0] = f_now;

			for (int j = 1; j <= max_order; ++j)
			{
				std::fill(phi[j].begin(), phi[j].end(), 0.0);
			}
		}

		std::vector<double> y_pred = y;
		for (int j = 0; j < order; ++j)
		{
			for (size_t i = 0; i < n; ++i)
			{
				y_pred[i] += h * G_Coeff[j] * phi[j][i];
			}
		}

		if (!Options.NonNegative.empty())
		{
			for (int idx : Options.NonNegative)
			{
				if (idx >= 0 && idx < static_cast<int>(n))
				{
					y_pred[idx] = std::max(y_pred[idx], 0.0);
				}
			}
		}

		std::vector<double> f_pred = F(t + h, y_pred);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		std::vector<std::vector<double>> d(max_order + 2, std::vector<double>(n));
		d[0] = f_pred;
		for (int j = 1; j <= order + 1; ++j)
		{
			for (size_t i = 0; i < n; ++i)
			{
				d[j][i] = d[j - 1][i] - phi[j - 1][i];
			}
		}

		std::vector<double> y_corr = y_pred;
		for (size_t i = 0; i < n; ++i)
		{
			y_corr[i] += h * G_Coeff[order] * d[order][i];
		}

		if (!Options.NonNegative.empty())
		{
			for (int idx : Options.NonNegative)
			{
				if (idx >= 0 && idx < static_cast<int>(n))
				{
					y_corr[idx] = std::max(y_corr[idx], 0.0);
				}
			}
		}

		std::vector<double> err_vec(n);
		for (size_t i = 0; i < n; ++i)
		{
			err_vec[i] = h * G_Coeff[order + 1] * d[order + 1][i];
		}

		double err = 0.0;
		for (size_t i = 0; i < n; ++i)
		{
			if (!std::isfinite(y_corr[i]))
			{
				throw std::runtime_error("ode113: non-finite solution at t = " + std::to_string(t + h));
			}

			double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y[i]), std::abs(y_corr[i])));
			err = std::max(err, std::abs(err_vec[i]) / tol);
		}

		if (err <= 1.0 || std::abs(h) < 16.0 * eps * std::abs(t))
		{
			if (Options.Stats)
			{
				res.stats.nsteps++;
			}

			std::vector<double> y_prev = y;
			std::vector<double> f_prev = f_now;

			f_now = F(t + h, y_corr);
			if (Options.Stats)
			{
				res.stats.nfevals++;
			}

			t += h;
			y = y_corr;

			std::vector<std::vector<double>> phi_old = phi;
			phi[0] = f_now;

			for (int j = 1; j <= max_order + 1; ++j)
			{
				for (size_t i = 0; i < n; ++i)
				{
					phi[j][i] = phi[j - 1][i] - phi_old[j - 1][i];
				}
			}

			if (return_specific_points)
			{
				while (tspan_idx < tspan.size() && ((dir > 0.0 && tspan[tspan_idx] <= t) || (dir < 0.0 && tspan[tspan_idx] >= t)))
				{
					double s = (tspan[tspan_idx] - (t - h)) / h;
					s = std::clamp(s, 0.0, 1.0);
					res.Time.push_back(tspan[tspan_idx]);
					res.Result.push_back(hermite_interp(s, h, y_prev, f_prev, y, f_now));
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
						res.Time.push_back(t - h + s * h);
						res.Result.push_back(hermite_interp(s, h, y_prev, f_prev, y, f_now));
					}
				}

				res.Time.push_back(t);
				res.Result.push_back(y);
			}

			if (Options.Events)
			{
				auto ev_prev = (*Options.Events)(t - h, y_prev);
				auto ev_next = (*Options.Events)(t, y);

				for (size_t ev_idx = 0; ev_idx < ev_prev.value.size(); ++ev_idx)
				{
					double v0 = ev_prev.value[ev_idx];
					double v1 = ev_next.value[ev_idx];

					bool crossed = false;
					if (v0 * v1 < 0.0)
					{
						int direction = (ev_idx < ev_prev.direction.size())	? ev_prev.direction[ev_idx] : 0;
						if (direction == 0)
						{
							crossed = true;
						}
						else if (direction > 0 && v1 > v0)
						{
							crossed = true;
						}
						else if (direction < 0 && v1 < v0)
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
							auto y_mid = hermite_interp(s_mid, h, y_prev, f_prev, y, f_now);
							auto ev_mid = (*Options.Events)(t - h + s_mid * h, y_mid);
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
						auto y_ev = hermite_interp(s_ev, h, y_prev, f_prev, y, f_now);
						double t_ev = t - h + s_ev * h;

						res.Te.push_back(t_ev);
						res.Ye.push_back(y_ev);
						res.Ie.push_back(static_cast<int>(ev_idx) + 1);

						bool terminal = (ev_idx < ev_prev.isterminal.size()) ? ev_prev.isterminal[ev_idx] : false;
						if (terminal)
						{
							return res;
						}
					}
				}
			}

			double h_opt = h * std::pow(0.5 / std::max(err, 1e-10), 1.0 / (order + 2));

			if (order < max_order)
			{
				std::vector<double> err_up_vec(n);
				for (size_t i = 0; i < n; ++i)
				{
					err_up_vec[i] = h * G_Coeff[order + 2] * d[order + 2][i];
				}
				double err_up = 0.0;
				for (size_t i = 0; i < n; ++i)
				{
					double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y_prev[i]), std::abs(y[i])));
					err_up = std::max(err_up, std::abs(err_up_vec[i]) / tol);
				}

				double h_opt_up = h * std::pow(0.5 / std::max(err_up, 1e-10), 1.0 / (order + 3));

				if (h_opt_up > h_opt)
				{
					h_opt = h_opt_up;
					order++;
				}
			}

			if (order > 1)
			{
				std::vector<double> err_down_vec(n);
				for (size_t i = 0; i < n; ++i)
				{
					err_down_vec[i] = h * G_Coeff[order] * d[order][i];
				}
				double err_down = 0.0;
				for (size_t i = 0; i < n; ++i)
				{
					double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y_prev[i]), std::abs(y[i])));
					err_down = std::max(err_down, std::abs(err_down_vec[i]) / tol);
				}

				double h_opt_down= h * std::pow(0.5 / std::max(err_down, 1e-10), 1.0 / (order + 1));

				if (h_opt_down > h_opt)
				{
					h_opt = h_opt_down;
					order--;
				}
			}

			h_old = h;
			h = dir * std::min(std::abs(h_opt), h_max);
			h = dir * std::max(std::abs(h), 16.0 * eps * std::abs(t));

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
			done = false;

			h *= 0.5;
			order = 1;
			phi[0] = f_now;
			for (int j = 1; j <= max_order; ++j)
			{
				std::fill(phi[j].begin(), phi[j].end(), 0.0);
			}
		}
	}

	return res;
}