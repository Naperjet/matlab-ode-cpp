#include "ode-suite.hpp"
#include <string>

namespace
{
	constexpr long double G[15] =
	{
		1.0L,
		1.0L / 2.0L,
		5.0L / 12.0L,
		3.0L / 8.0L,
		251.0L / 720.0L,
		95.0L / 288.0L,
		19087.0L / 60480.0L,
		5257.0L / 17280.0L,
		1070017.0L / 3628800.0L,
		25713.0L / 89600.0L,
		26842253.0L / 95800320.0L,
		4777223.0L / 17418240.0L,
		703604254357.0L / 2615348736000.0L,
		106364763817.0L / 402361344000.0L,
		1166309819657.0L / 4483454976000.0L
	};

	std::vector<long double> hermite_dense(long double s, long double h, const std::vector<long double>& y0, const std::vector<long double>& f0, const std::vector<long double>& y1, const std::vector<long double>& f1)
	{
		long double s2 = s * s;
		long double s3 = s2 * s;
		long double h00 = 2.0L * s3 - 3.0L * s2 + 1.0L;
		long double h10 = s3 - 2.0L * s2 + s;
		long double h01 = -2.0L * s3 + 3.0L * s2;
		long double h11 = s3 - s2;

		size_t n = y0.size();
		std::vector<long double> interp(n);
		for (size_t i = 0; i < n; ++i)
			interp[i] = h00 * y0[i] + h01 * y1[i] + h * (h10 * f0[i] + h11 * f1[i]);
		return interp;
	}

	long double estimate_initial_step(long double t0, long double tf, const std::vector<long double>& y0, const std::vector<long double>& f0, long double RelTol, long double AbsTol)
	{
		size_t n = y0.size();
		long double ny = 0.0L;
		long double nf = 0.0L;

		for (size_t i = 0; i < n; ++i)
		{
			long double sc = std::max(AbsTol, RelTol * std::abs(y0[i]));
			ny = std::max(ny, std::abs(y0[i]) / sc);
			nf = std::max(nf, std::abs(f0[i]) / sc);
		}

		long double h0 = (ny < 1e-5L || nf < 1e-5L) ? 1e-6L : 0.01L * (ny / nf);
		return std::min(h0, std::abs(tf - t0));
	}
}

ODE_Res_Long ode113_long(const ODE_Func_Long& F, const std::vector<long double>& tspan, const std::vector<long double>& Init_Value, const ODE_Set_Long& Options)
{
	if (tspan.size() < 2)
	{
		throw std::invalid_argument("ode113(long ver): tspan must contain at least two points.");
	}

	bool increasing = (tspan[1] >= tspan[0]);
	for (size_t i = 1; i < tspan.size(); ++i)
	{
		if (increasing && tspan[i] < tspan[i - 1])
		{
			throw std::invalid_argument("ode113(long ver): tspan must be strictly increasing");
		}
		if (!increasing && tspan[i] > tspan[i - 1])
		{
			throw std::invalid_argument("ode113(long ver): tspan must be strictly decreasing");
		}
	}

	ODE_Res_Long res;
	const bool return_specific_points = (tspan.size() > 2);
	const long double t0 = tspan.front();
	const long double tf = tspan.back();
	long double t = t0;

	std::vector<long double> y = Init_Value;
	size_t n = y.size();

	if (n == 0)
	{
		throw std::invalid_argument("ode113(long ver): initial state cannot be empty");
	}

	res.Time.push_back(t);
	res.Result.push_back(y);

	long double dir = (tf >= t0) ? 1.0L : -1.0L;

	long double h;
	if (Options.Init_Step > 0.0L)
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

	long double h_max = (Options.Max_Step > 0.0L) ? Options.Max_Step : std::abs(tf - t0) / 10.0L;
	h_max = std::abs(h_max);
	const long double eps = std::numeric_limits<long double>::epsilon();
	const size_t max_iter = 1'000'000;

	const int max_order = 13;

	std::vector<std::vector<long double>> phi(max_order + 2, std::vector<long double>(n, 0.0L));

	std::vector<long double> f_now = F(t, y);
	if (Options.Stats)
	{
		res.stats.nfevals++;
	}
	phi[0] = f_now;

	int order = 1;
	long double h_old = h;
	size_t tspan_idx = 1;
	size_t iter = 0;

	bool completed = false;

	while (!completed)
	{
		if (++iter > max_iter)
		{
			throw std::runtime_error("ode113(long ver): maximum iterations reached");
		}

		if ((dir > 0.0L && t + h > tf) || (dir < 0.0L && t + h < tf))
		{
			h = tf - t;
			completed = true;
		}

		long double h_ratio = std::abs(h / h_old);
		if ((h_ratio < 0.5L || h_ratio > 2.0L) && order > 1)
		{
			order = 1;
			phi[0] = f_now;
			for (int j = 1; j <= max_order + 1; ++j)
				std::fill(phi[j].begin(), phi[j].end(), 0.0L);
		}

		std::vector<long double> y_pred = y;
		for (int j = 0; j < order; ++j)
		{
			for (size_t i = 0; i < n; ++i)
			{
				y_pred[i] += h * G[j] * phi[j][i];
			}
		}

		if (!Options.NonNegative.empty())
		{
			for (int idx : Options.NonNegative)
			{
				if (idx >= 0 && idx < static_cast<int>(n))
				{
					y_pred[idx] = std::max(y_pred[idx], 0.0L);
				}
			}
		}

		std::vector<long double> f_pred = F(t + h, y_pred);
		if (Options.Stats) 
		{
			res.stats.nfevals++;
		}

		std::vector<std::vector<long double>> d(max_order + 2, std::vector<long double>(n));
		d[0] = f_pred;
		for (int j = 1; j <= order + 1; ++j)
		{
			for (size_t i = 0; i < n; ++i)
			{
				d[j][i] = d[j - 1][i] - phi[j - 1][i];
			}
		}

		std::vector<long double> y_corr = y_pred;
		for (size_t i = 0; i < n; ++i)
		{
			y_corr[i] += h * G[order] * d[order][i];
		}

		if (!Options.NonNegative.empty())
		{
			for (int idx : Options.NonNegative)
			{
				if (idx >= 0 && idx < static_cast<int>(n))
				{
					y_corr[idx] = std::max(y_corr[idx], 0.0L);
				}
			}
		}

		std::vector<long double> err_vec(n);
		for (size_t i = 0; i < n; ++i)
			err_vec[i] = h * G[order + 1] * d[order + 1][i];

		long double err = 0.0L;
		for (size_t i = 0; i < n; ++i)
		{
			if (!std::isfinite(y_corr[i]))
			{
				throw std::runtime_error("ode113(long ver): non-finite solution at t = " + std::to_string(static_cast<double>(t + h)));
			}

			long double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y[i]), std::abs(y_corr[i])));
			err = std::max(err, std::abs(err_vec[i]) / tol);
		}

		if (err <= 1.0L || std::abs(h) <= 16.0L * eps * std::abs(t))
		{
			if (Options.Stats)
			{
				res.stats.nsteps++;
			}

			std::vector<long double> y_prev = y;
			std::vector<long double> f_prev = f_now;

			f_now = F(t + h, y_corr);
			if (Options.Stats)
			{
				res.stats.nfevals++;
			}

			t += h;
			y = y_corr;

			std::vector<std::vector<long double>> phi_old = phi;
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
				while (tspan_idx < tspan.size() && ((dir > 0.0L && tspan[tspan_idx] <= t) || (dir < 0.0L && tspan[tspan_idx] >= t)))
				{
					long double s = (tspan[tspan_idx] - (t - h)) / h;
					s = std::clamp(s, 0.0L, 1.0L);
					res.Time.push_back(tspan[tspan_idx]);
					res.Result.push_back(hermite_dense(s, h, y_prev, f_prev, y, f_now));
					tspan_idx++;
				}
			}
			else
			{
				if (Options.Refine > 1)
				{
					for (int j = 1; j < Options.Refine; ++j)
					{
						long double s = static_cast<long double>(j) / static_cast<long double>(Options.Refine);
						res.Time.push_back(t - h + s * h);
						res.Result.push_back(hermite_dense(s, h, y_prev, f_prev, y, f_now));
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
					long double v0 = ev_prev.value[ev_idx];
					long double v1 = ev_next.value[ev_idx];

					bool crossed = false;
					if (v0 * v1 < 0.0L)
					{
						int direction = (ev_idx < ev_prev.direction.size()) ? ev_prev.direction[ev_idx] : 0;
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
						long double s_lo = 0.0L, s_hi = 1.0L;
						long double v_lo = v0, v_hi = v1;

						for (int b = 0; b < 30; ++b)
						{
							long double s_mid = 0.5L * (s_lo + s_hi);
							auto y_mid = hermite_dense(s_mid, h, y_prev, f_prev, y, f_now);
							auto ev_mid = (*Options.Events)(t - h + s_mid * h, y_mid);
							long double v_mid = ev_mid.value[ev_idx];

							if (v_lo * v_mid <= 0.0L)
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

						long double s_ev = 0.5L * (s_lo + s_hi);
						auto y_ev = hermite_dense(s_ev, h, y_prev, f_prev, y, f_now);
						long double t_ev = t - h + s_ev * h;

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

			long double h_opt = h * std::pow(0.5L / std::max(err, 1e-18L), 1.0L / (order + 2));
			int next_order = order;

			if (order < max_order)
			{
				std::vector<long double> err_up_vec(n);
				for (size_t i = 0; i < n; ++i)
				{
					err_up_vec[i] = h * G[order + 2] * d[order + 2][i];
				}

				long double err_up = 0.0L;
				for (size_t i = 0; i < n; ++i)
				{
					long double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y_prev[i]), std::abs(y[i])));
					err_up = std::max(err_up, std::abs(err_up_vec[i]) / tol);
				}

				long double h_opt_up = h * std::pow(0.5L / std::max(err_up, 1e-18L), 1.0L / (order + 3));
				if (h_opt_up > h_opt)
				{
					h_opt = h_opt_up;
					next_order = order + 1;
				}
			}

			if (order > 1)
			{
				std::vector<long double> err_down_vec(n);
				for (size_t i = 0; i < n; ++i)
				{
					err_down_vec[i] = h * G[order] * d[order][i];
				}

				long double err_down = 0.0L;
				for (size_t i = 0; i < n; ++i)
				{
					long double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y_prev[i]), std::abs(y[i])));
					err_down = std::max(err_down, std::abs(err_down_vec[i]) / tol);
				}

				long double h_opt_down = h * std::pow(0.5L / std::max(err_down, 1e-18L), 1.0L / (order + 1));
				if (h_opt_down > h_opt)
				{
					h_opt = h_opt_down;
					next_order = order - 1;
				}
			}

			if (next_order > 9 && err < 1e-13L)
			{
				next_order = 9;
			}

			if (next_order > 7 && err < 1e-15L)
			{
				next_order = 7;
			}

			h_old = h;
			h = dir * std::min(std::abs(h_opt), h_max);
			h = dir * std::max(std::abs(h), 16.0L * eps * std::abs(t));
			order = next_order;

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
			completed = false;

			h *= 0.5L;
			order = 1;
			phi[0] = f_now;
			for (int j = 1; j <= max_order + 1; ++j)
			{
				std::fill(phi[j].begin(), phi[j].end(), 0.0L);
			}
		}
	}

	return res;
}