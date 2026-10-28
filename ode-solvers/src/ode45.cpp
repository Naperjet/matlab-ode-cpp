#include "ode-suite.hpp"

namespace
{
	std::vector<double> dopri5_dense(double theta, double h, const std::vector<double>& y, const std::vector<double>& y_next, const std::vector<std::vector<double>>& k)
	{

		static const double P[7][4] =
		{
			{  1.0,                     -8048581381.0 / 2820520608.0,   8663915743.0 / 2820520608.0,  -12715105075.0 / 11282082432.0 },
			{  0.0,                      0.0,                          0.0,                          0.0                       },
			{  0.0,                    131558114200.0 / 32700410799.0, -68118460800.0 / 10900136933.0,  87487479700.0 / 32700410799.0 },
			{  0.0,                     -1754552775.0 / 470086768.0,    14199869525.0 / 1410260304.0,   -10690763975.0 / 1880347072.0 },
			{  0.0,                    127303824393.0 / 49829197408.0, -318862633887.0 / 49829197408.0, 701980252875.0 / 199316789632.0 },
			{  0.0,                     -282668133.0 / 205662961.0,      2019193451.0 / 616988883.0,    -1453857185.0 / 822651844.0   },
			{  0.0,                     40617522.0 / 29380423.0,        -110615467.0 / 29380423.0,       69997945.0 / 29380423.0      }
		};

		const double th2 = theta * theta;
		const double th3 = th2 * theta;
		const double th4 = th3 * theta;

		double coeff[7];
		for (int s = 0; s < 7; ++s)
		{
			coeff[s] = P[s][0] * theta + P[s][1] * th2 + P[s][2] * th3 + P[s][3] * th4;
		}

		size_t n = y.size();
		std::vector<double> interp(n);
		for (size_t i = 0; i < n; ++i)
		{
			double acc = 0.0;
			for (int s = 0; s < 7; ++s)
			{
				acc += coeff[s] * k[s][i];
			}
			interp[i] = y[i] + h * acc;
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

ODE_Res ode45(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options)
{
	if (tspan.size() < 2)
	{
		throw std::invalid_argument("ode45: tspan must contain at least two elements");
	}

	bool increasing = (tspan[1] >= tspan[0]);
	for (size_t i = 1; i < tspan.size(); ++i)
	{
		if (increasing && tspan[i] < tspan[i - 1])
		{
			throw std::invalid_argument("ode45: Tspan is not strictly increasing.");
		}

		if (!increasing && tspan[i] > tspan[i - 1])
		{
			throw std::invalid_argument("ode45: Tspan is not strictly increasing.");
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
		throw std::invalid_argument("ode45: Initial state is empty");
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

	double h_max = (Options.Max_Step > 0.0) ? Options.Max_Step : std::abs(tf - t0) / 10;
	h_max = std::abs(h_max);

	const double eps = std::numeric_limits<double>::epsilon();
	const size_t max_iter = 1'000'000;

	std::vector<std::vector<double>> k(7, std::vector<double>(n));
	std::vector<double> y_next(n);

	k[0] = F(t, y);
	if (Options.Stats)
	{
		res.stats.nfevals++;
	}

	size_t tspan_idx = 1;
	size_t iter = 0;
	double err_old = 1e-4;

	while ((direction > 0.0 && t < tf) || (direction < 0.0 && t > tf))
	{
		if (++iter > max_iter)
		{
			throw std::runtime_error("ode45: Consider stiff solvers. Maximum iteration count exceeded");
		}

		if ((direction > 0.0 && t + h > tf) || (direction < 0.0 && t + h < tf))
		{
			h = tf - t;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * (1.0 / 5.0) * k[0][i];
		}
		k[1] = F(t + h * (1.0 / 5.0), y_next);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * ((3.0 / 40.0) * k[0][i] + (9.0 / 40.0) * k[1][i]);
		}
		k[2] = F(t + h * (3.0 / 10.0), y_next);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * ((44.0 / 45.0) * k[0][i] - (56.0 / 15.0) * k[1][i] + (32.0 / 9.0) * k[2][i]);
		}
		k[3] = F(t + h * (4.0 / 5.0), y_next);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * ((19372.0 / 6561.0) * k[0][i] - (25360.0 / 2187.0) * k[1][i] + (64448.0 / 6561.0) * k[2][i] - (212.0 / 729.0) * k[3][i]);
		}
		k[4] = F(t + h * (8.0 / 9.0), y_next);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * ((9017.0 / 3168.0) * k[0][i] - (355.0 / 33.0) * k[1][i] + (46732.0 / 5247.0) * k[2][i] + (49.0 / 176.0) * k[3][i] - (5103.0 / 18656.0) * k[4][i]);
		}
		k[5] = F(t + h, y_next);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y_next[i] = y[i] + h * ((35.0 / 384.0) * k[0][i] + (500.0 / 1113.0) * k[2][i] + (125.0 / 192.0) * k[3][i] - (2187.0 / 6784.0) * k[4][i] + (11.0 / 84.0) * k[5][i]);
		}
		k[6] = F(t + h, y_next);
		if (Options.Stats)
		{
			res.stats.nfevals++;
		}

		double max_err_ratio = 0.0;
		for (size_t i = 0; i < n; i++)
		{
			double error = h * ((35.0 / 384.0 - 5179.0 / 57600.0) * k[0][i] + (500.0 / 1113.0 - 7571.0 / 16695.0) * k[2][i] + (125.0 / 192.0 - 393.0 / 640.0) * k[3][i] - (2187.0 / 6784.0 - 92097.0 / 339200.0) * k[4][i] + (11.0 / 84.0 - 187.0 / 2100.0) * k[5][i] - (1.0 / 40.0) * k[6][i]);

			double tol = std::max(Options.Abs_Tol, Options.Rel_Tol * std::max(std::abs(y[i]), std::abs(y_next[i])));

			max_err_ratio = std::max(max_err_ratio, std::abs(error) / tol);
		}

		bool stop_early = false;

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
					double theta = (tspan[tspan_idx] - t) / h;
					res.Time.push_back(tspan[tspan_idx]);
					res.Result.push_back(dopri5_dense(theta, h, y, y_next, k));
					++tspan_idx;
				}
			}
			else
			{
				if (Options.Refine > 1)
				{
					for (int j = 1; j < Options.Refine; ++j)
					{
						double theta = static_cast<double>(j) / static_cast<double>(Options.Refine);

						res.Time.push_back(t + theta * h);
						res.Result.push_back(dopri5_dense(theta, h, y, y_next, k));
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
						int dir = (ev_idx < ev_prev.direction.size())
							? ev_prev.direction[ev_idx] : 0;
						if (dir == 0) crossed = true;
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
						double th_lo = 0.0, th_hi = 1.0;
						double v_lo = v0, v_hi = v1;

						for (int b = 0; b < 30; ++b)
						{
							double th_mid = 0.5 * (th_lo + th_hi);
							auto y_mid = dopri5_dense(th_mid, h, y, y_next, k);
							auto ev_mid = (*Options.Events)(t + th_mid * h, y_mid);
							double v_mid = ev_mid.value[ev_idx];

							if (v_lo * v_mid <= 0.0)
							{
								th_hi = th_mid;
								v_hi = v_mid;
							}
							else
							{
								th_lo = th_mid;
								v_lo = v_mid;
							}
						}

						double th_ev = 0.5 * (th_lo + th_hi);
						auto y_ev = dopri5_dense(th_ev, h, y, y_next, k);
						double t_ev = t + th_ev * h;

						res.Te.push_back(t_ev);
						res.Ye.push_back(y_ev);
						res.Ie.push_back(static_cast<int>(ev_idx) + 1);

						bool terminal = (ev_idx < ev_prev.isterminal.size())
							? ev_prev.isterminal[ev_idx] : false;
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
			k[0] = k[6];

			err_old = std::max(max_err_ratio, 1e-10);

			if (return_specific_points && tspan_idx >= tspan.size())
			{
				stop_early = true;
			}
		}
		else
		{
			if (Options.Stats)
			{
				res.stats.nfailed++;
			}
		}

		double h_scale;
		if (max_err_ratio < 1e-15)
		{
			h_scale = 5.0;
		}
		else
		{
			h_scale = 0.8 * std::pow(1.0 / std::max(max_err_ratio, 1e-10), 0.2) * std::pow(err_old / std::max(max_err_ratio, 1e-10), 0.08);
		}

		h_scale = std::clamp(h_scale, 0.1, 5.0);
		h = direction * std::max(std::abs(h * h_scale), 16.0 * eps * std::abs(t));
		h = direction * std::min(std::abs(h), h_max);

		if (std::abs(h) < 16.0 * eps * std::abs(t) && max_err_ratio > 1.0)
		{
			throw std::runtime_error("ode45: step size below minimum. Consider a stiff solver.");
		}

		if (stop_early)
		{
			break;
		}
	}

	return res;
}