#include "Lin-Alg.hpp"

#include<algorithm>
#include<limits>

mat Jacob_make(RHS_Fcn f, double t, const vec& y, const vec& f_val)
{
	size_t n = y.size();
	mat Jacob(static_cast<int>(n), static_cast<int>(n));

	const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

	vec y_p = y;
	double* y_ptr = y_p.data.data();
	const double* y_orig_ptr = y.data.data();
	double* jacob_ptr = Jacob.data.data();
	const double* f_val_ptr = f_val.data.data();

	for (size_t j = 0; j < n; ++j)
	{
		double yj = y_orig_ptr[j];
		double step = eps * std::max(std::abs(yj), 1.0);
		y_ptr[j] = yj + step;

		vec f_p = f(t, y_p);
		const double* f_p_ptr = f_p.data.data();

		double inv_step = 1.0 / step;
		size_t col_offset = j * n;

		for (size_t i = 0; i < n; ++i)
		{
			jacob_ptr[col_offset + i] = (f_p_ptr[i] - f_val_ptr[i]) * inv_step;
		}

		y_ptr[j] = yj;
	}

	return Jacob;
}

mat Jacob_make(const std::function<std::vector<double>(double, const std::vector<double>&)>& f, double t, const std::vector<double>& y, const std::vector<double>& f_val)
{
	RHS_Fcn f_vec = [&f](double t_in, const vec& y_in) -> vec
		{
			return vec(f(t_in, y_in.data));
		};

	vec y_vec(y);
	vec f_val_vec(f_val);

	return(Jacob_make(f_vec, t, y_vec, f_val_vec));
}