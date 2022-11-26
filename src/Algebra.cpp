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

    std::vector<double> difference(const std::vector<double> &a, const std::vector<double> &b)
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

    void summ(std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b)
    {
        size_type n = a.size();

        for (size_type i = 0; i < n; ++i)
        {
            for (size_type j = 0; j < n; ++j)
            {
                a.at(i).at(j) += b.at(i).at(j);
            }
        }
    }

    void difference(std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b)
    {
        for (size_type i = 0; i < a.size(); ++i)
        {
            for (size_type j = 0; j < a.at(0).size(); ++j)
            {
                a.at(i).at(j) -= b.at(i).at(j);
            }
        }
    }

    void multiply(std::vector<std::vector<double>> &a, double value)
    {
        for (auto& i: a)
        {
            for (auto& j: i)
            {
                j *= value;
            }
        }
    }

    std::vector<std::vector<double>>
    multiply(const std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b)
    {
        size_type n = a.size();
        std::vector<std::vector<double>> res = zero(n);

        // TODO optimize this
        for (size_type i = 0; i < n; ++i)
        {
            for (size_type j = 0; j < n; ++j)
            {
                for (size_type k = 0; k < n; ++k)
                {
                    res.at(i).at(j) = a.at(i).at(k) * b.at(k).at(j);
                }
            }
        }

        return res;
    }

    std::vector<std::vector<double>> one(size_type n)
    {
        std::vector<std::vector<double>> res(n);

        for (auto& i: res)
        {
            i = std::vector<double>(n);
        }

        for (size_type i = 0; i < n; ++i)
        {
            res[i][i] = 1.0;
        }

        return res;
    }

    std::vector<std::vector<double>> zero(size_type n)
    {
        std::vector<std::vector<double>> res(n);

        for (auto& i: res)
        {
            i = std::vector<double>(n);
        }

        return res;
    }

    std::vector<double> summ(const std::vector<double> &a, const std::vector<double> &b)
    {
        size_type n = a.size();
        std::vector<double> result(n);

        for (size_type i = 0; i < n; ++i)
        {
            result[i] = a[i] + b[i];
        }

        return result;
    }

    std::vector<double> multiply(const std::vector<double> &a, double value)
    {
        size_type n = a.size();
        std::vector<double> result(n);

        for (size_type i = 0; i < n; ++i)
        {
            result[i] = a[i] * value;
        }

        return result;
    }
}
