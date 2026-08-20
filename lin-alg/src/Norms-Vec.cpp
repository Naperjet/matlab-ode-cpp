#include "Lin-Alg.hpp"

double dot(const vec& v1, const vec& v2)
{
	assert(v1.size() == v2.size() && "Dot Product Dimension Mismatch");
	double sum = 0.0;
	size_t n = v1.size();

	for (size_t i = 0; i < n; ++i)
	{
		sum += v1[i] * v2[i];
	}

	return sum;
}

double norm_inf(const vec& v)
{
	double max_val = 0.0;
	size_t n = v.size();

	for (size_t i = 0; i < n; ++i)
	{
		double abs_val = std::abs(v[i]);
		if (abs_val > max_val)
		{
			max_val = abs_val;
		}
	}

	return max_val;
}

double norm_wts(const vec& v, const vec& weights)
{
	assert(v.size() == weights.size() && "Norm Dimension Mismatch");
	double sum = 0.0;
	size_t n = v.size();

	for (size_t i = 0; i < n; ++i)
	{
		double prod = v[i] * weights[i];
		sum += prod * prod;
	}

	return std::sqrt(sum / n);
}

vec operator*(double scalar, const vec& v)
{
	return scalar * v;
}