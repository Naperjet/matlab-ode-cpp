#include "Lin-Alg.hpp"

mat mat::operator+(const mat& other) const
{
	assert(rows == other.rows && cols == other.cols);
	mat result(rows, cols);

	size_t n = data.size();

	const double* a = data.data();
	const double* b = other.data.data();

	double* r = result.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		r[i] = a[i] + b[i];
	}

	return result;
}

mat mat::operator-(const mat& other) const
{
	assert(rows == other.rows && cols == other.cols);
	mat result(rows, cols);

	size_t n = data.size();

	const double* a = data.data();
	const double* b = other.data.data();

	double* r = result.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		r[i] = a[i] - b[i];
	}

	return result;
}

mat mat::operator*(double scalar) const
{
	mat result(rows, cols);
	size_t n = data.size();
	const double* a = data.data();
	double* r = result.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		r[i] = a[i] * scalar;
	}

	return result;
}

mat& mat::operator+=(const mat& other)
{
	assert(rows == other.rows && cols == other.cols);
	size_t n = data.size();
	double* a = data.data();
	const double* b = other.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		a[i] += b[i];
	}

	return *this;
}

mat& mat::operator-=(const mat& other)
{
	assert(rows == other.rows && cols == other.cols);
	size_t n = data.size();
	double* a = data.data();
	const double* b = other.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		a[i] -= b[i];
	}

	return *this;
}

mat& mat::operator*=(double scalar)
{
	size_t n = data.size();
	double* a = data.data();

	for (size_t i = 0; i < n; ++i)
	{
		a[i] *= scalar;
	}

	return *this;
}
