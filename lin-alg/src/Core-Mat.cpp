#include "Lin-Alg.hpp"

mat::mat(int r, int c, double initial_val) : data(static_cast<size_t>(r) * c, initial_val), rows(r), cols(c) {}

double& mat::operator()(int r, int c)
{
	assert(r >= 0 && r < rows && "Row index out of bounds");
	assert(c >= 0 && c < cols && "Column index out of bounds");

	return data[static_cast<size_t>(c) * rows + r];
}

const double& mat::operator()(int r, int c) const
{
	assert(r >= 0 && r < rows && "Row index out of bounds");
	assert(c >= 0 && c < cols && "Column index out of bounds");

	return data[static_cast<size_t>(c) * rows + r];
}