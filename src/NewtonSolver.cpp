//
// Created by Severin on 18.11.2022.
//

#include "NewtonSolver.h"

#include <cmath>
#include <utility>

NewtonSolver::NewtonSolver(const Context& context)
{
    m_context = context;
}

std::vector<double> NewtonSolver::solve_newton(std::function<std::vector<double>(std::vector<double>)> &f,
                                               const vector &initial_guess, uint32_t max_iterations)
{
    size_type n = initial_guess.size();

    vector previous = initial_guess;
    vector next = vector(n);

    for (uint32_t i = 0; i < max_iterations; ++i)
    {
        // n by n + 1 matrix: (derivatives | derivatives * previous - f(previous))
        matrix m;
        create(m, n, n + 1);

        // Find all derivatives (last column stay intact)
        derive(m, f, previous);

        // Find right part and put into m (compute last column)
        plug_vector(m, subtract(multiply(m, previous), f(previous)));

        solve(m, next);

        if (residual(previous, next) < m_context.newton_precision)
        {
            break;
        }

        std::swap(next, previous);
    }

    vector output = next;

    return output;
}


double NewtonSolver::residual(const vector &a, const vector &b)
{
    double res = 0;
    size_type n = a.size();

    for (size_type i = 0; i < n; ++i)
    {
        res += (a[i] - b[i]) * (a[i] - b[i]);
    }

    return res;
}

void NewtonSolver::derive(matrix &m, std::function<std::vector<double>(std::vector<double>)> &f, const vector &x) const
{
    size_type n = m.size();

    for (size_type column = 0; column < n; ++column)
    {
        vector l = x;
        vector r = x;
        l[column] += m_context.newton_derive_step;
        r[column] -= m_context.newton_derive_step;

        vector d = divide(subtract(f(l), f(r)), 2 * m_context.newton_derive_step);

        for (size_type row = 0; row < n; ++row)
        {
            m[row][column] = d[row];
        }
    }
}

void NewtonSolver::solve(matrix &mx, vector &res)
{
    size_type n = mx.size();
    vector_int x_order = vector_int(n);

    for (size_type i = 0; i < n; i++)
    {
        x_order.at(i) = i;
    }

    for (size_type j = 0; j < n; j++)
    {
        // Serch for biggest element and change order //

        size_type max_i = 0;
        size_type max_j = 0;
        double max = 0;

        for (size_type z = j; z < n; z++)
        {
            for (size_type t = j; t < n; t++)
            {
                if (fabs(mx.at(z).at(t)) > max)
                {
                    max = fabs(mx.at(z).at(t));
                    max_i = z;
                    max_j = t;
                }
            }
        }

        for (size_type i = 0; i < n; i++)
        {
            std::swap(mx.at(i).at(max_j), mx.at(i).at(j));
        }
        std::swap(mx.at(max_i), mx.at(j));
        std::swap(x_order.at(j), x_order.at(max_j));

        ////////////////////////////////////////////////

        for (size_type i = j + 1; i < n; i++)
        {
            double div = mx.at(i).at(j) / mx.at(j).at(j);

            for (size_type k = 0; k < n + 1; k++)
            {
                mx.at(i).at(k) = mx.at(i).at(k) - div * mx.at(j).at(k);
            }
        }
    }

    res.at(n - 1) = mx.at(n - 1).at(n) / mx.at(n - 1).at(n - 1);

    for (int i = (int) n - 2; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < (int) n; j++)
        {
            sum = sum + mx.at(i).at(j) * res.at(j);
        }
        res.at(i) = (mx.at(i).at(n) - sum) / mx.at(i).at(i);
    }

    // Restore the order of res //

    for (size_type i = 0; i < n; i++)
    {
        size_type next = i;

        while (x_order.at(next) != n)
        {

            std::swap(res.at(i), res.at(x_order.at(next)));

            size_type temp = x_order.at(next);

            x_order.at(next) = n;
            next = temp;
        }
    }
}

NewtonSolver::NewtonSolver()
{
    m_context = Context();
}
