#include "ode-suite.hpp"
#include "Lin-Alg.hpp"
#include <deque>

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

vec deval(double t_query, const ODE_Res& sol)
{
    const auto& t = sol.Time;
    const auto& r = sol.Result;
    size_t num_points = t.size();

    if (num_points < 4)
    {
        throw std::runtime_error("Deval Requires at least 4 data points");
    }

    if (t_query <= t.front())
    {
        return vec(r.front());
    }

    if (t_query >= t.back())
    {
        return vec(r.back());
    }

    auto it = std::lower_bound(t.begin(), t.end(), t_query);
    size_t idx = std::distance(t.begin(), it);

    size_t start_idx = (idx > 1) ? idx - 2 : 0;
    if (start_idx + 3 >= num_points)
    {
        start_idx = num_points - 4;
    }

    double t0 = t[start_idx];
    double t1 = t[start_idx + 1];
    double t2 = t[start_idx + 2];
    double t3 = t[start_idx + 3];

    vec y0(r[start_idx]);
    vec y1(r[start_idx + 1]);
    vec y2(r[start_idx + 2]);
    vec y3(r[start_idx + 3]);

    double L0 = ((t_query - t1) * (t_query - t2) * (t_query - t3)) / ((t0 - t1) * (t0 - t2) * (t0 - t3));
    double L1 = ((t_query - t0) * (t_query - t2) * (t_query - t3)) / ((t1 - t0) * (t1 - t2) * (t1 - t3));
    double L2 = ((t_query - t0) * (t_query - t1) * (t_query - t3)) / ((t2 - t0) * (t2 - t1) * (t2 - t3));
    double L3 = ((t_query - t0) * (t_query - t1) * (t_query - t2)) / ((t3 - t0) * (t3 - t1) * (t3 - t2));

    return (y0 * L0) + (y1 * L1) + (y2 * L2) + (y3 * L3);
}

ODE_Res ode15s_core(ODEFunc F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, ODE_Set Options)
{
    ODE_Res res;
    double t0 = tspan.front();
    double tf = tspan.back();

    size_t n_eq = Init_Value.size();
    vec y(Init_Value);

    res.Time.reserve(10000);
    res.Result.reserve(10000);

    res.Time.push_back(t0);
    res.Result.push_back(y.data);

    RHS_Fcn f_wrap = [&] (double t, const vec& y_in) -> vec
        {
            return vec(F(t, y_in.data));
        };
    
    double t = t0;
    double h = (Options.Init_Step > 0) ? Options.Init_Step : 1e-4;
    double h_max = (Options.Max_Step > 0) ? Options.Max_Step : (tf - t0) / 10;

    size_t current_order = 1;
    const size_t max_order = 5;
    const int max_newton_iter = 7;

    std::deque<vec> history;
    history.push_front(y);

    mat J(n_eq, n_eq);
    mat A(n_eq, n_eq);
    std::vector<size_t> pivots;

    bool requires_jacobian_update = true;
    bool requires_matrix_factorization = true;

    while (t < tf)
    {
        if (t + h > tf)
        {
            h = tf - t;
        }

        bool step_accepted = false;
        vec y_next = y;

        while (!step_accepted)
        {
            if (requires_jacobian_update)
            {
                vec f_val = f_wrap (t, y);
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

                catch(...)
                {
                    h *= 0.2;
                    requires_matrix_factorization = true;
                    continue;
                }

                requires_matrix_factorization = false;
            }

            vec y_iter = y_next;
            bool newton_converged = false;
            vec delta(n_eq);

            vec history_sum(n_eq, 0.0);
            const auto& alphas = BDF_Alpha[current_order - 1];
            for (size_t i = 0; i < current_order; i++)
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

                if (norm_inf(delta) < 0.1 * Options.Rel_Tol)
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

                current_order = 1;
                history.resize(1);
                continue;
            }

            vec error_weights(n_eq);
            for (size_t i = 0; i < n_eq; i++)
            {
                error_weights[i] = 1.0 / (Options.Abs_Tol + Options.Rel_Tol * std::max(std::abs(y[i]), std::abs(y_next[i])));
            }

            double error_norm = norm_wts(delta, error_weights);

            if (error_norm <= 1.0)
            {
                t += h;
                y = y_next;
                res.Time.push_back(t);
                res.Result.push_back(y.data);
                step_accepted = true;

                history.push_front(y);
                if (history.size() > max_order)
                {
                    history.pop_back();
                }

                if (current_order < max_order && history.size() == current_order + 1)
                {
                    current_order++;
                    requires_matrix_factorization = true;
                }

                double scale = 0.9 * std::pow(1.0 / std::max(error_norm, 1e-10), 1.0 / (current_order + 1));
                scale = std::min(std::max(scale, 0.2), 2.0);

                double h_new = std::min(h * scale, h_max);
                if (std::abs(h_new - h) > 0.1 * h)
                {
                    h = h_new;
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
                    history.pop_back();
                }
            }
        }
    }

    return res;
}

ODE_Res ode15s(ODEFunc F, const std::vector<double>& tspan, const std::vector<double>& Init_Value, ODESet Options)
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