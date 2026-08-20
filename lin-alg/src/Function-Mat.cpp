#include "Lin-Alg.hpp"

vec mat_vec_multiply(const mat& m, const vec& v)
{
	size_t cols = static_cast<size_t>(m.cols);
	size_t rows = static_cast<size_t>(m.rows);
	assert(cols == v.size() && "Matrix and Vector Multiplication Dimension Mismatch");

	vec result(rows, 0.0);

	double* res_ptr = result.data.data();
	const double* mat_ptr = m.data.data();
	const double* vec_ptr = v.data.data();

	for (size_t c = 0; c < cols; ++c)
	{
		double v_val = vec_ptr[c];
		size_t col_offset = c * rows;

		for (size_t r = 0; r < rows; ++r)
		{
			res_ptr[r] += mat_ptr[col_offset + r] * v_val;
		}
	}

	return result;
}

mat identity(int size)
{
	mat I(size, size, 0.0);

	for (int i = 0; i < size; ++i)
	{
		I(i, i) = 1.0;
	}

	return I;
}

mat operator*(double scalar, const mat& m)
{
	return m * scalar;
}