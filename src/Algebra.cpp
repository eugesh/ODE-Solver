//
// Created by Severin on 18.11.2022.
//

#include "Algebra.h"

namespace algebra
{
    std::vector<double> multiply(const std::vector<std::vector<double>> &m, const std::vector<double> &x)
    {
        size_type n = m.size();
        std::vector<double> result(n);

        for (size_type i = 0; i < n; ++i)
        {
            for (size_type j = 0; j < n; ++j)
            {
                result[i] += m[i][j] * x[j];
            }
        }

        return result;
    }

    std::vector<double> subtract(const std::vector<double> &a, const std::vector<double> &b)
    {
        size_type n = a.size();
        std::vector<double> result(n);

        for (size_type i = 0; i < n; ++i)
        {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    std::vector<double> divide(const std::vector<double> &a, double value)
    {
        size_type n = a.size();
        std::vector<double> result(n);

        for (size_type i = 0; i < n; ++i)
        {
            result[i] = a[i] / value;
        }

        return result;
    }

    void plug_vector(std::vector<std::vector<double>> &m, const std::vector<double> &x)
    {
        size_type n = m.size();

        for (size_type i = 0; i < n; ++i)
        {
            m[i][n] = x[i];
        }
    }
}
