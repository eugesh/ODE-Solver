//
// Created by Severin on 24.11.2022.
//

#include "Utils.h"

Utils::Utils(const Context &context)
{
    m_context = context;
}

void Utils::derive(std::vector<std::vector<double>> &mx, std::function<std::vector<double>(std::vector<double>)> &f,
                   const std::vector<double> &x) const
{
    algebra::size_type n = mx.size();

    std::vector<double> l = x;
    std::vector<double> r = x;

    for (algebra::size_type column = 0; column < n; ++column)
    {
        l[column] += m_context.derive_step;
        r[column] -= m_context.derive_step;

        std::vector<double> d = algebra::divide(algebra::difference(f(l), f(r)), 2 * m_context.derive_step);

        for (algebra::size_type row = 0; row < n; ++row)
        {
            mx[row][column] = d[row];
        }

        l[column] -= m_context.derive_step;
        r[column] += m_context.derive_step;
    }
}

void Utils::solve(std::vector<std::vector<double>> &mx, std::vector<double> &res)
{
    algebra::size_type n = mx.size();
    std::vector<algebra::size_type> x_order(n);

    for (algebra::size_type i = 0; i < n; i++)
    {
        x_order[i] = i;
    }

    for (algebra::size_type j = 0; j < n; j++)
    {
        // Search for biggest element and change order //

        algebra::size_type max_i = 0;
        algebra::size_type max_j = 0;
        double max = 0;

        for (algebra::size_type z = j; z < n; z++)
        {
            for (algebra::size_type t = j; t < n; t++)
            {
                if (fabs(mx[z][t]) > max)
                {
                    max = fabs(mx[z][t]);
                    max_i = z;
                    max_j = t;
                }
            }
        }

        for (algebra::size_type i = 0; i < n; i++)
        {
            std::swap(mx[i][max_j], mx[i][j]);
        }
        std::swap(mx[max_i], mx[j]);
        std::swap(x_order[j], x_order[max_j]);

        ////////////////////////////////////////////////

        for (algebra::size_type i = j + 1; i < n; i++)
        {
            double div = mx[i][j] / mx[j][j];

            for (algebra::size_type k = 0; k < n + 1; k++)
            {
                mx[i][k] = mx[i][k] - div * mx[j][k];
            }
        }
    }

    res[n - 1] = mx[n - 1][n] / mx[n - 1][n - 1];

    for (int i = (int) n - 2; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < (int) n; j++)
        {
            sum = sum + mx[i][j] * res[j];
        }
        res[i] = (mx[i][n] - sum) / mx[i][i];
    }

    // Restore the order of res //

    for (algebra::size_type i = 0; i < n; i++)
    {
        algebra::size_type next = i;

        while (x_order[next] != n)
        {

            std::swap(res[i], res[x_order[next]]);

            algebra::size_type temp = x_order[next];

            x_order[next] = n;
            next = temp;
        }
    }
}

double Utils::residual(const std::vector<double> &a, const std::vector<double> &b)
{
    double res = 0;
    algebra::size_type n = a.size();

    for (algebra::size_type i = 0; i < n; ++i)
    {
        res += (a[i] - b[i]) * (a[i] - b[i]);
    }

    return sqrt(res);
}
