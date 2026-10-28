#include "ode-suite.hpp"
#include "Lin-Alg.hpp"
#include <deque>

namespace
{
	const std::vector<std::vector<double>> BDF_Alpha =
	{
		{1.0},
		{4.0 / 3.0, -1.0 / 3.0},
		{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0},
		{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0},
		{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}
	};

	const std::vector<double> BDF_Beta =
	{
		1.0, 2.0 / 3.0, 6.0 / 11.0, 12.0 / 25.0, 60.0 / 137.0
	};

	double estimate_initial_step(double t0, double tf, const vec& y0, const vec& f0, double reltol, double abstol)
	{
		size_t n = y0.data.size();
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

	vec newton_divided_diff_eval(double t_query, const std::vector<double>& times, const std::vector<vec>& values)
	{
		size_t k = times.size();
		if (k == 0)
		{
			throw std::runtime_error("newton_divided_diff_eval: no interpolation nodes supplied");
		}

		size_t n = values[0].data.size();
		std::vector<vec> coef = values;

		for (size_t j = 1; j < k; ++j)
		{
			for (size_t i = k; i-- > j; )
			{
				for (size_t c = 0; c < n; ++c)
				{
					coef[i][c] = (coef[i][c] - coef[i - 1][c]) / (times[i] - times[i - j]);
				}
			}
		}

		vec result = coef[k - 1];
		for (size_t idx = k - 1; idx-- > 0; )
		{
			for (size_t c = 0; c < n; ++c)
			{
				result[c] = coef[idx][c] + (t_query - times[idx]) * result[c];
			}
		}

		return result;
	}

	double wrms_norm(const vec& v, const vec& y, const vec& y_next, double abstol, double reltol)
	{
		size_t n = y.data.size();
		double sum = 0.0;

		for (size_t i = 0; i < n; ++i)
		{
			double sc = abstol + reltol * std::max(std::abs(y[i]), std::abs(y_next[i]));
			double vi = v[i] / sc;
			sum += vi * vi;
		}

		return std::sqrt(sum / n);
	}
}

vec deval(double t_query, const ODE_Res& sol)
{
	const auto& t = sol.Time;
	const auto& r = sol.Result;
	size_t num_points = t.size();

	if (num_points < 2)
	{
		throw std::runtime_error("deval requires at least 2 data points");
	}

	if (t_query <= t.front())
	{
		return vec(r.front());
	}

	if (t_query >= t.back())
	{
		return vec(r.back());
	}

	size_t window = std::min<size_t>(4, num_points);

	auto it = std::lower_bound(t.begin(), t.end(), t_query);
	size_t idx = std::distance(t.begin(), it);

	size_t start_idx;
	if (window == 4)
	{
		start_idx = (idx > 1) ? idx - 2 : 0;
		if (start_idx + 3 >= num_points)
		{
			start_idx = num_points - 4;
		}
	}
	else
	{
		start_idx = 0;
		if (start_idx + (window - 1) >= num_points)
		{
			start_idx = num_points - window;
		}
	}

	std::vector<double> nodes_t(window);
	std::vector<vec> nodes_y(window);
	for (size_t i = 0; i < window; ++i)
	{
		nodes_t[i] = t[start_idx + i];
		nodes_y[i] = vec(r[start_idx + i]);
	}

	return newton_divided_diff_eval(t_query, nodes_t, nodes_y);
}

ODE_Res ode15s_core(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options)
{
	if (tspan.size() < 2)
	{
		throw std::invalid_argument("ode15s: tspan must have at least two elements");
	}

	bool increasing = (tspan[1] > tspan[0]);
	for (size_t i = 1; i < tspan.size(); ++i)
	{
		if (increasing && tspan[i] < tspan[i - 1])
		{
			throw std::invalid_argument("ode15s: tspan must be strictly increasing");
		}
		if (!increasing && tspan[i] > tspan[i - 1])
		{
			throw std::invalid_argument("ode15s: tspan must be strictly decreasing");
		}
	}

	double t0 = tspan.front();
	double tf = tspan.back();
	double direction = (tf >= t0) ? 1.0 : -1.0;

	size_t n_eq = Init_Value.size();
	if (n_eq == 0)
	{
		throw std::invalid_argument("ode15s: initial state cannot be empty");
	}

	ODE_Res res;
	vec y(Init_Value);

	res.Time.reserve(10000);
	res.Result.reserve(10000);

	res.Time.push_back(t0);
	res.Result.push_back(y.data);

	RHS_Fcn f_wrap = [&](double t_in, const vec& y_in) -> vec
		{
			return vec(F(t_in, y_in.data));
		};

	double t = t0;

	double h;
	if (Options.Init_Step > 0)
	{
		h = Options.Init_Step;
	}
	else
	{
		vec f0 = f_wrap(t, y);
		h = estimate_initial_step(t0, tf, y, f0, Options.Rel_Tol, Options.Abs_Tol);
	}
	h = direction * std::abs(h);

	double h_max = (Options.Max_Step > 0) ? Options.Max_Step : std::abs(tf - t0) / 10;
	h_max = std::abs(h_max);

	const double eps = std::numeric_limits<double>::epsilon();
	const double h_min = 16.0 * eps * std::max({ std::abs(t0), std::abs(tf), 1.0 });

	size_t current_order = 1;
	const size_t max_order = 5;
	const int max_newton_iter = 7;

	std::deque<vec> history;
	std::deque<double> time_history;
	history.push_front(y);
	time_history.push_front(t);

	mat J(n_eq, n_eq);
	mat A(n_eq, n_eq);
	std::vector<size_t> pivots;

	bool requires_jacobian_update = true;
	bool requires_matrix_factorization = true;

	while ((direction > 0.0 && t < tf) || (direction < 0.0 && t > tf))
	{
		if ((direction > 0.0 && t + h > tf) || (direction < 0.0 && t + h < tf))
		{
			h = tf - t;
		}

		bool step_accepted = false;
		vec y_next = y;

		while (!step_accepted)
		{
			if (std::abs(h) < h_min)
			{
				throw std::runtime_error("ode15s: step size too small; failed to converge to acceptable tolerances");
			}

			if (requires_jacobian_update)
			{
				vec f_val = f_wrap(t, y);
				J = Jacob_make(f_wrap, t, y, f_val);
				requires_jacobian_update = false;
				requires_matrix_factorization = true;
			}

			if (requires_matrix_factorization)
			{
				double beta = BDF_Beta[current_order - 1];
				A = identity(n_eq);
				for (size_t c = 0; c < n_eq; c++)
				{
					for (size_t r = 0; r < n_eq; r++)
					{
						A(r, c) = A(r, c) - h * beta * J(r, c);
					}
				}

				try
				{
					lu_factor(A, pivots);
				}

				catch (...)
				{
					h *= 0.2;
					requires_jacobian_update = true;
					requires_matrix_factorization = true;
					continue;
				}

				requires_matrix_factorization = false;
			}

			const auto& alphas = BDF_Alpha[current_order - 1];

			size_t corr_avail = std::min(current_order, history.size());

			size_t pred_avail = std::min(current_order + 1, history.size());

			std::vector<double> nodes_t(pred_avail);
			std::vector<vec> nodes_y(pred_avail);
			for (size_t i = 0; i < pred_avail; ++i)
			{
				nodes_t[i] = time_history[i];
				nodes_y[i] = history[i];
			}
			vec y_iter = newton_divided_diff_eval(t + h, nodes_t, nodes_y);
			const vec y_predictor = y_iter;

			bool newton_converged = false;
			vec delta(n_eq);

			vec history_sum(n_eq, 0.0);
			for (size_t i = 0; i < corr_avail; ++i)
			{
				history_sum += (history[i] * alphas[i]);
			}

			double h_beta = h * BDF_Beta[current_order - 1];
			vec G(n_eq);

			for (int iter = 0; iter < max_newton_iter; iter++)
			{
				vec f_next = f_wrap(t + h, y_iter);

				double* G_ptr = G.data.data();
				const double* y_iter_ptr = y_iter.data.data();
				const double* hist_ptr = history_sum.data.data();
				const double* f_next_ptr = f_next.data.data();

				for (size_t i = 0; i < n_eq; i++)
				{
					G_ptr[i] = y_iter_ptr[i] - hist_ptr[i] - (f_next_ptr[i] * h_beta);
				}

				delta = lu_solve(A, pivots, G);
				y_iter -= delta;

				double delta_wrms = wrms_norm(delta, y, y_iter, Options.Abs_Tol, Options.Rel_Tol);

				if (delta_wrms < 0.1)
				{
					newton_converged = true;
					y_next = y_iter;
					break;
				}
			}

			if (!newton_converged)
			{
				h *= 0.5;
				requires_jacobian_update = true;
				requires_matrix_factorization = true;

				if (current_order > 1)
				{
					--current_order;
				}
				continue;
			}

			vec lte = y_next - y_predictor;
			double error_norm = wrms_norm(lte, y, y_next, Options.Abs_Tol, Options.Rel_Tol);

			if (error_norm <= 1.0)
			{
				size_t order_used = current_order;

				t += h;
				y = y_next;
				res.Time.push_back(t);
				res.Result.push_back(y.data);
				step_accepted = true;

				history.push_front(y);
				time_history.push_front(t);
				if (history.size() > max_order + 1)
				{
					history.pop_back();
					time_history.pop_back();
				}

				if (current_order < max_order && history.size() >= current_order + 2 && error_norm < 0.5)
				{
					current_order++;
					requires_matrix_factorization = true;
				}

				double scale = 0.9 * std::pow(1.0 / std::max(error_norm, 1e-10), 1.0 / (order_used + 1));
				scale = std::min(std::max(scale, 0.2), 2.0);

				double h_new = std::min(std::abs(h * scale), h_max);
				if (std::abs(h_new - std::abs(h)) > 0.1 * std::abs(h))
				{
					h = direction * h_new;
					requires_matrix_factorization = true;
				}
			}

			else
			{
				double scale = 0.9 * std::pow(1.0 / error_norm, 1.0 / (current_order + 1));
				scale = std::max(scale, 0.1);
				h *= scale;
				requires_matrix_factorization = true;

				if (current_order > 1)
				{
					current_order--;
				}
			}
		}
	}

	return res;
}

ODE_Res ode15s(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options)
{
	ODE_Res raw_res = ode15s_core(F, tspan, Init_Value, Options);

	if (tspan.size() <= 2)
	{
		return raw_res;
	}

	ODE_Res dense_res;
	dense_res.Time = tspan;
	dense_res.Result.reserve(tspan.size());

	for (double target_time : tspan)
	{
		vec interp = deval(target_time, raw_res);
		dense_res.Result.push_back(interp.data);
	}

	return dense_res;
}