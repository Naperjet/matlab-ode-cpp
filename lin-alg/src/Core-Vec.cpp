#include "Lin-Alg.hpp"

vec::vec(size_t size, double initial_val) :data(size, initial_val) {}

vec::vec(const std::vector<double>& init_data) :data(init_data) {}

size_t vec::size() const
{
    return data.size();
}

double& vec::operator[](int i)
{
    return data[i];
}

const double& vec::operator[](int i) const
{
    return data[i];
}