#include "Lin-Alg.hpp"

vec vec::operator+(const vec& other) const
{
	size_t n = size();
	assert(n == other.size() && "Vector Addition Dimension Mismatch");
	vec result(n);

	double* res_ptr = result.data.data();
	const double* this_ptr = data.data();
	const double* other_ptr = other.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		res_ptr[i] = this_ptr[i] + other_ptr[i];
	}

	return result;
}

vec vec::operator-(const vec& other) const
{
	size_t n = size();
	assert(n == other.size() && "Vector Subtraction Dimension Mismatch");
	vec result(n);

	double* res_ptr = result.data.data();
	const double* this_ptr = data.data();
	const double* other_ptr = other.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		res_ptr[i] = this_ptr[i] - other_ptr[i];
	}

	return result;
}

vec vec::operator*(double scalar) const
{
	size_t n = size();
	vec result(n);

	double* res_ptr = result.data.data();
	const double* this_ptr = data.data();

	for (size_t i = 0; i < n; ++i)
	{
		res_ptr[i] = this_ptr[i] * scalar;
	}

	return result;
}

vec& vec::operator+=(const vec& other)
{
	size_t n = size();
	assert(n == other.size() && "Vector += Dimension Mismatch");

	double* this_ptr = data.data();
	const double* other_ptr = other.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		this_ptr[i] += other_ptr[i];
	}

	return *this;
}

vec& vec::operator-=(const vec& other)
{
	size_t n = size();
	assert(n == other.size() && "Vector -= Dimension Mismatch");

	double* this_ptr = data.data();
	const double* other_ptr = other.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		this_ptr[i] -= other_ptr[i];
	}

	return *this;
}

vec& vec::operator*=(double scalar)
{
	size_t n = size();

	double* this_ptr = data.data();

	for (size_t i = 0; i < n; ++i)
	{
		this_ptr[i] *= scalar;
	}

	return *this;
}