#include "Lin-Alg.hpp"
#include<algorithm>

void lu_factor(mat& A, std::vector<size_t>& pivot_indices)
{
	size_t n = A.rows;
	assert(n == static_cast<size_t>(A.cols) && "LU factorization needs a square matrix");

	pivot_indices.resize(n);
	for (size_t i = 0; i < n; ++i)
	{
		pivot_indices[i] = i;
	}

	double* A_ptr = A.data.data();

	for (size_t k = 0; k < n; ++k)
	{
		double max_val = 0.0;
		size_t pivot_row = k;

		for (size_t i = k; i < n; ++i)
		{
			double abs_val = std::abs(A_ptr[k * n + i]);
			if (abs_val > max_val)
			{
				max_val = abs_val;
				pivot_row = i;
			}
		}

		const double threshold = std::numeric_limits<double>::epsilon() * n;
		if (max_val < threshold)
		{
			throw std::runtime_error("Singular Matrix Detected");
		}

		if (pivot_row != k)
		{
			std::swap(pivot_indices[k], pivot_indices[pivot_row]);

			for (size_t j = 0; j < n; ++j)
			{
				std::swap(A_ptr[j * n + k], A_ptr[j * n + pivot_row]);
			}
		}

		double A_kk = A_ptr[k * n + k];
		size_t col_k_offset = k * n;

		for (size_t i = k + 1; i < n; ++i)
		{
			A_ptr[col_k_offset + i] /= A_kk;
		}

		for (size_t j = k + 1; j < n; ++j)
		{
			double A_kj = A_ptr[j * n + k];
			size_t col_j_offset = j * n;

			for (size_t i = k + 1; i < n; ++i)
			{
				A_ptr[col_j_offset + i] -= A_ptr[col_k_offset + i] * A_kj;
			}
		}
	}
}

vec lu_solver(const mat& LU, const std::vector<size_t>& pivot_indices, const vec& b)
{
	size_t n = LU.rows;
	vec x(n);

	const double* LU_ptr = LU.data.data();
	double* x_ptr = x.data.data();
	const double* b_ptr = b.data.data();

	for (size_t i = 0; i < n; ++i)
	{
		x_ptr[i] = b_ptr[pivot_indices[i]];
	}

	for (size_t j = 0; j < n; ++j)
	{
		double x_j = x_ptr[j];
		size_t col_j_offset = j * n;

		for (size_t i = j + 1; i < n; ++i)
		{
			x_ptr[i] -= LU_ptr[col_j_offset + i] * x_j;
		}
	}

	for (size_t j = n; j-- > 0;)
	{
		size_t col_j_offset = j * n;
		x_ptr[j] /= LU_ptr[col_j_offset + j];
		double x_j = x_ptr[j];

		for (size_t i = 0; i < j; ++i)
		{
			x_ptr[i] -= LU_ptr[col_j_offset + i] * x_j;
		}
	}

	return x;
}