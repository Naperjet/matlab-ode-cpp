#include "ode-suite.hpp"
#include "Lin-Alg.hpp"

namespace
{
	const double d = (1.0 / (2.0 + std::sqrt(2.0)));
	const double e32 = (6.0 + std::sqrt(2.0));

	std::vector<double> hermite_dense(double theta, double h, const std::vector<double>& y, const std::vector<double>& y_next, const std::vector<double>& f, const std::vector<double>& f_next)
	{
		double th2 = theta * theta;
		double th3 = th2 * theta;

		double h0 = 2.0 * th3 - 3.0 * th2 + 1.0;
		double h1 = -2.0 * th3 + 3.0 * th2;
		double h2 = th3 - 2.0 * th2 + theta;
		double h3 = th3 - th2;

		size_t n = y.size();
		std::vector<double> interp(n);

		for (size_t i = 0; i < n; ++i)
		{
			interp[i] = h0 * y[i] + h1 * y_next[i] + h2 * h * f[i] + h3 * h * f_next[i];
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

	double wrms_norm(const std::vector<double>& v, const std::vector<double>& y, const std::vector<double>& y_next, double abstol, double reltol)
	{
		size_t n = v.size();
		double sum = 0.0;

		for (size_t i = 0; i < n; ++i)
		{
			double sc = abstol + reltol * std::max(std::abs(y[i]), std::abs(y_next[i]));
			double vi = v[i] / sc;
			sum += vi * vi;
		}

		return std::sqrt(sum / n);
	}

	void enforce_nonneg(const std::vector<int>& nonneg_idx, std::vector<double>& y, size_t n)
	{
		for (int idx : nonneg_idx)
		{
			if (idx >= 0 && static_cast<size_t>(idx) < n && y[idx] < 0.0)
			{
				y[idx] = 0.0;
			}
		}
	}
}

ODE_Res ode23s(const ODE_Func& F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, const ODE_Set& Options)
{
	if (tspan.size() < 2)
	{
		throw std::invalid_argument("ode23s: tspan must contain at least two elements");
	}

	bool increasing = (tspan[1] > tspan[0]);
	for (size_t i = 1; i < tspan.size(); ++i)
	{
		if (increasing && tspan[i] < tspan[i - 1])
		{
			throw std::invalid_argument("ode23s: tspan must be strictly increasing");
		}

		if (!increasing && tspan[i] > tspan[i - 1])
		{
			throw std::invalid_argument("ode23s: tspan must be strictly decreasing");
		}
	}

	const bool return_specific_points = (tspan.size() > 2);
	const double t0 = tspan.front();
	const double tf = tspan.back();

	double direction = (tf >= t0) ? 1.0 : -1.0;

	const size_t n = Init_Value.size();
	if (n == 0)
	{
		throw std::invalid_argument("ode23s: Initial state cannot be empty");
	}


	ODE_Res res;
	res.Time.push_back(t0);
	res.Result.push_back(Init_Value);

	double t = t0;
	std::vector<double> y = Init_Value;

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

	const double h_min = 16.0 * eps * std::max({ std::abs(t0), std::abs(tf), 1.0 });

	const size_t max_iter = 1'000'000;

	int n_int = static_cast<int>(n);
	mat J(n_int, n_int);
	mat W(n_int, n_int);
	std::vector<size_t> pivots;

	bool needs_jac = true;

	bool needs_fac = true;

	std::vector<double> f1(n);
	std::vector<double> f2(n);

	std::vector<double> k1(n);
	std::vector<double> k2(n);
	std::vector<double> k3(n);

	std::vector<double> G2(n);
	std::vector<double> G3(n);

	std::vector<double> y_tmp(n);
	std::vector<double> y_next(n);
	std::vector<double> delta(n);

	RHS_Fcn f_wrap = [&](double t_in, const vec& y_in) -> vec
		{
			if (Options.Stats)
			{
				res.stats.nfevals++;
			}
			return vec(F(t_in, y_in.data));
		};

	std::vector<double> f_current = F(t, y);
	if (Options.Stats)
	{
		res.stats.nfevals++;
	}

    size_t tspan_idx = 1;
    size_t iter = 0;
    size_t steps_since_jac = 0;        
    const size_t max_steps_per_jac = 5; 
    double h_at_last_jac = std::abs(h); 

    while ((direction > 0.0 && t < tf) || (direction < 0.0 && t > tf))
    {
        if (++iter > max_iter)
        {
            throw std::runtime_error("ode23s: Maximum iteration exceeded");
        }

        if ((direction > 0.0 && t + h > tf) || (direction < 0.0 && t + h < tf))
        {
            h = tf - t;
        }

        bool step_accepted = false;

        while (!step_accepted)
        {
            if (std::abs(h) < h_min)
            {
                throw std::runtime_error("ode23s: step size too small; failed to converge to acceptable tolerances");
            }

            bool h_changed_drastically = (std::abs(h) > 2.0 * h_at_last_jac) ||
                (std::abs(h) < 0.5 * h_at_last_jac);

            if (needs_jac || steps_since_jac >= max_steps_per_jac || h_changed_drastically)
            {
                vec y_vec(y);
                vec f0_vec(f_current);
                J = Jacob_make(f_wrap, t, y_vec, f0_vec);
                needs_jac = false;
                needs_fac = true;
                steps_since_jac = 0;
                h_at_last_jac = std::abs(h);
            }

            if (needs_fac)
            {
                W = identity(n_int);
                double hd = h * d;

                for (size_t c = 0; c < n; ++c)
                {
                    for (size_t r = 0; r < n; ++r)
                    {
                        W(static_cast<int>(r), static_cast<int>(c)) -= hd * J(static_cast<int>(r), static_cast<int>(c));
                    }
                }

                try
                {
                    lu_factor(W, pivots);
                }
                catch (...)
                {
                    h *= 0.2;
                    if (Options.Stats)
                    {
                        res.stats.nfailed++;
                    }
                    needs_fac = true;
                    needs_jac = true;
                    continue;
                }

                needs_fac = false;
            }

            vec f0_vec = f_current;
            vec k1_vec = lu_solve(W, pivots, f0_vec);

            for (size_t i = 0; i < n; ++i)
            {
                k1[i] = k1_vec[i];
            }

            double h05 = 0.5 * h;
            for (size_t i = 0; i < n; ++i)
            {
                y_tmp[i] = y[i] + h05 * k1[i];
            }

            f1 = F(t + h05, y_tmp);
            if (Options.Stats)
            {
                res.stats.nfevals++;
            }

            for (size_t i = 0; i < n; ++i)
            {
                G2[i] = f1[i] - k1[i] + f_current[i];
            }

            vec G2_vec(G2);
            vec k2_vec = lu_solve(W, pivots, G2_vec);

            for (size_t i = 0; i < n; ++i)
            {
                k2[i] = k2_vec[i];
            }

            for (size_t i = 0; i < n; ++i)
            {
                y_next[i] = y[i] + h * k2[i];
            }

            f2 = F(t + h, y_next);
            if (Options.Stats)
            {
                res.stats.nfevals++;
            }

            for (size_t i = 0; i < n; ++i)
            {
                G3[i] = f2[i] - e32 * k2[i] - 2.0 * (k1[i] + k2[i]) + f_current[i];
            }

            vec G3_vec(G3);
            vec k3_vec = lu_solve(W, pivots, G3_vec);

            for (size_t i = 0; i < n; ++i)
            {
                k3[i] = k3_vec[i];
            }

            for (size_t i = 0; i < n; ++i)
            {
                delta[i] = h * k3[i];
            }

            double error_norm = wrms_norm(delta, y, y_next, Options.Abs_Tol, Options.Rel_Tol);

            if (error_norm <= 1.0)
            {
                if (Options.Stats)
                {
                    res.stats.nsteps++;
                }

                steps_since_jac++;

                double t_prev = t;
                std::vector<double> y_prev = y;
                std::vector<double> f_prev = f_current;

                t += h;
                y = std::move(y_next);
                f_current = std::move(f2);

                if (!Options.NonNegative.empty())
                {
                    std::vector<double> y_before_clamp = y;
                    enforce_nonneg(Options.NonNegative, y, n);
                    if (y != y_before_clamp)
                    {
                        f_current = F(t, y);
                        if (Options.Stats)
                        {
                            res.stats.nfevals++;
                        }
                    }
                }

                if (return_specific_points)
                {
                    while (tspan_idx < tspan.size() && ((direction > 0.0 && tspan[tspan_idx] <= t) || (direction < 0.0 && tspan[tspan_idx] >= t)))
                    {
                        double theta = (tspan[tspan_idx] - t_prev) / h;
                        theta = std::clamp(theta, 0.0, 1.0); 
                        res.Time.push_back(tspan[tspan_idx]);
                        res.Result.push_back(hermite_dense(theta, h, y_prev, y, f_prev, f_current));
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
                            res.Time.push_back(t_prev + theta * h);
                            res.Result.push_back(hermite_dense(theta, h, y_prev, y, f_prev, f_current));
                        }
                    }

                    res.Time.push_back(t);
                    res.Result.push_back(y);
                }

                if (Options.Events)
                {
                    auto ev_prev = (*Options.Events)(t_prev, y_prev);
                    auto ev_next = (*Options.Events)(t, y);

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
                            double th_lo = 0.0, th_hi = 1.0;
                            double v_lo = v0, v_hi = v1;

                            for (int b = 0; b < 30; ++b)
                            {
                                double th_mid = 0.5 * (th_lo + th_hi);
                                auto y_mid = hermite_dense(th_mid, h, y_prev, y, f_prev, f_current);
                                auto ev_mid = (*Options.Events)(t_prev + th_mid * h, y_mid);
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
                            auto y_ev = hermite_dense(th_ev, h, y_prev, y, f_prev, f_current);
                            double t_ev = t_prev + th_ev * h;

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

                double scale = 0.9 * std::pow(1.0 / std::max(error_norm, 1e-10), 1.0 / 3.0);
                scale = std::clamp(scale, 0.2, 2.0);

                double h_new = std::min(std::abs(h * scale), h_max);
                if (std::abs(h_new - std::abs(h)) > 0.1 * std::abs(h))
                {
                    h = direction * h_new;
                    needs_fac = true;
                }

                step_accepted = true;

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

                double scale = 0.9 * std::pow(1.0 / error_norm, 1.0 / 3.0);
                scale = std::max(scale, 0.1);
                h *= scale;
                needs_fac = true;

                if (scale < 0.3)
                {
                    needs_jac = true;
                }
            }
        }
    }

    return res;
}