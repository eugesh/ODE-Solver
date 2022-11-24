//
// Created by Severin on 18.11.2022.
//

#include "NewtonSolver.h"

#include <cmath>

NewtonSolver::NewtonSolver(const Context &context)
{
    m_context = context;
    m_utils = Utils(context);

    // Create matrix of derivatives
    size_type n = context.x_0.size();

    m = std::vector<std::vector<double>>(n);
    for (auto &item: m)
    {
        item = std::vector<double>(n + 1);
    }
}

std::vector<double> NewtonSolver::solve_newton(std::function<std::vector<double>(std::vector<double>)> &f,
                                               const std::vector<double> &initial_guess, uint32_t max_iterations)
{
    size_type n = initial_guess.size();

    std::vector<double> tmp = initial_guess;
    std::vector<double> output(n);

    for (uint32_t i = 0; i < max_iterations; ++i)
    {
        // Find all derivatives (last column stay intact)
        m_utils.derive(m, f, tmp);

        // Find right part and put into m (compute last column)
        plug_vector(m, subtract(multiply(m, tmp), f(tmp)));

        // Solve system of linear equations
        m_utils.solve(m, output);

        if (m_utils.residual(tmp, output) < m_context.newton_precision)
        {
            break;
        }

        std::swap(output, tmp);
    }

    return output;
}

