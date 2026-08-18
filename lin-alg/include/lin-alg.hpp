#pragma once

#include <vector>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <functional>


struct vec
{
	std::vector<double> data;

	vec() = default;

	vec(size_t size, double initial_val = 0.0);
	vec(const std::vector<double>& init_data);

	size_t size() const;

	double& operator[](size_t i);
	const double& operator[](size_t i) const;

	vec operator+(const vec& other) const;
	vec operator-(const vec& other) const;
	vec operator*(double scalar) const;

	vec& operator+=(const vec& other);
	vec& operator-=(const vec& other);
	vec& operator*=(double scalar);
};

struct mat
{
	std::vector<double> data;

	mat() = default;

	int rows;
	int cols;

	mat(int r, int c, double initial_val = 0.0);

	double& operator() (int r, int c);
	const double& operator() (int r, int c) const;

	mat operator+(const mat& other) const;
	mat operator-(const mat& other) const;
	mat operator*(double scalar) const;

	mat& operator+=(const mat& other);
	mat& operator-=(const mat& other);
	mat& operator*=(double scalar);
};

double dot(const vec& v1, const vec& v2);
double norm_inf(const vec& v);

double norm_wts(const vec& v, const vec& weights);

vec mat_vec_multiply(const mat& m, const vec& v);
mat identity(int size);

void lu_factor(mat& A, std::vector<size_t>& pivot_indices);

vec lu_solve(const mat& LU, const std::vector<size_t>& pivot_indices, const vec& b);

mat Jacob_make(const std::function<std::vector<double>(double, const std::vector<double>&)>& f, double t, const std::vector<double>& y, const std::vector<double>& f_val);

using RHS_Fcn = std::function<vec(double, const vec&)>;

mat Jacob_make(RHS_Fcn f, double t, const vec& y, const vec& f_val);

mat operator*(double scalar, const mat& m);

vec operator*(double scalar, const vec& v);