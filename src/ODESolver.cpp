//
// Created by Severin on 08.11.2022.
//

#include "ODESolver.h"

#include <cmath>
#include <functional>
#include <utility>
#include <iostream>

void ODESolver::rk()
{
    if (!m_context.f)
    {
        std::cerr << "ODE is not defined!" << std::endl;
        std::exit(1);
    }

    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    while (t < m_context.t_end)
    {

        // TODO optimize this
        std::vector<double> k_1 = k1(t, x);
        std::vector<double> k_2 = k2(t, x, k_1);
        std::vector<double> k_3 = k3(t, x, k_2);
        std::vector<double> k_4 = k4(t, x, k_3);

        for (algebra::size_type i = 0; i < dim; ++i)
        {
            x.at(i) += 1.0 / 6.0 * (k_1.at(i) + 2 * k_2.at(i) + 2 * k_3.at(i) + k_4.at(i));
        }

        Result.emplace_back(t, x);

        t += m_context.h;
    }
}

std::vector<double> ODESolver::k1(double t, std::vector<double> x) const
{
    std::vector<double> res = m_context.f(t, std::move(x));

    for (auto &item: res)
    {
        item *= m_context.h;
    }

    return res;
}

std::vector<double> ODESolver::k2(double t, std::vector<double> x, std::vector<double> k_1) const
{
    for (algebra::size_type i = 0; i < x.size(); ++i)
    {
        x.at(i) += k_1.at(i) / 2;
    }

    std::vector<double> res = m_context.f(t + m_context.h / 2, std::move(x));

    for (auto &item: res)
    {
        item *= m_context.h;
    }

    return res;

}

std::vector<double> ODESolver::k3(double t, std::vector<double> x, std::vector<double> k_2) const
{
    for (algebra::size_type i = 0; i < x.size(); ++i)
    {
        x.at(i) += k_2.at(i) / 2;
    }

    std::vector<double> res = m_context.f(t + m_context.h / 2, std::move(x));

    for (auto &item: res)
    {
        item *= m_context.h;
    }

    return res;

}

std::vector<double> ODESolver::k4(double t, std::vector<double> x, std::vector<double> k_3) const
{
    for (algebra::size_type i = 0; i < x.size(); ++i)
    {
        x.at(i) += k_3.at(i);
    }

    std::vector<double> res = m_context.f(t + m_context.h, std::move(x));

    for (auto &item: res)
    {
        item *= m_context.h;
    }

    return res;

}

void ODESolver::ae()
{
    if (!m_context.f)
    {
        std::cerr << "ODE is not defined!" << std::endl;
        std::exit(1);
    }

    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    computeA();

    // Compute n initial values of x
    ComputeInitialAdamsX(t, x);

    // Run Adams method

    // Create n initial values of f
    std::vector<std::vector<double>> fs = ComputeInitialAdamsF();

    while (t < m_context.t_end)
    {

        for (algebra::size_type j = 0; j < m_context.adams_order; ++j)
        {
            for (algebra::size_type i = 0; i < dim; ++i)
            {
                x.at(i) += A.at(m_context.adams_order - 1 - j) * fs.at(j).at(i) * m_context.h;
            }
        }

        Result.emplace_back(t, x);

        t += m_context.h;

        for (uint32_t i = 0; i < fs.size() - 1; ++i)
        {
            fs.at(i) = std::move(fs.at(i + 1));
        }

        fs.at(fs.size() - 1) = m_context.f(t, x);
    }
}

void ODESolver::ai()
{
    if (!m_context.f)
    {
        std::cerr << "ODE is not defined!" << std::endl;
        std::exit(1);
    }

    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    computeB();

    // Compute n initial values of x
    ComputeInitialAdamsX(t, x);

    // Compute initial f values
    std::vector<std::vector<double>> fs = ComputeInitialAdamsF();

    // Run Adams method

    std::function<std::vector<double>(std::vector<double>)> f =
            [&fs, &t, this, &x](std::vector<double> y) -> std::vector<double>
            {
                std::vector<double> res = m_context.f(t, y);

                for (double &re: res)
                {
                    re *= B[0];
                }

                for (algebra::size_type j = 0; j < m_context.adams_order - 1; ++j)
                {
                    for (algebra::size_type i = 0; i < dim; ++i)
                    {
                        res[i] += B[m_context.adams_order - 1 - j] * fs[j + 1][i];
                    }
                }

                for (double &re: res)
                {
                    re *= m_context.h;
                }

                for (algebra::size_type i = 0; i < dim; ++i)
                {
                    res[i] += x[i] - y[i];
                }

                return res;
            };

    while (t < m_context.t_end)
    {
        x = m_newton_solver.solve_newton(f, x);
        Result.emplace_back(t, x);

        t += m_context.h;

        for (algebra::size_type i = 0; i < fs.size() - 1; ++i)
        {
            fs.at(i) = std::move(fs.at(i + 1));
        }

        fs.at(fs.size() - 1) = m_context.f(t, x);
    }
}

void ODESolver::computeA()
{
    A = std::vector<double>(m_context.adams_order);

    //TODO optimize this
    for (int32_t j = 0; j < (int32_t) m_context.adams_order; ++j)
    {
        A.at(j) = pow(-1.0, j) / (factorial(j) * factorial((int32_t) m_context.adams_order - 1 - j));

        A.at(j) *= integrate(integrandForA, j);
    }
}

void ODESolver::computeB()
{
    B = std::vector<double>(m_context.adams_order);

    //TODO optimize this
    for (int32_t j = -1; j < (int32_t) m_context.adams_order - 1; ++j)
    {
        B.at(j + 1) = pow(-1.0, j + 1) / (factorial(j + 1) * factorial((int32_t) m_context.adams_order - 2 - j));

        B.at(j + 1) *= integrate(integrandForB, j);
    }
}

int32_t ODESolver::factorial(int32_t x)
{
    int32_t res = 1;

    for (int32_t i = 2; i <= x; ++i)
    {
        res *= i;
    }

    return res;
}

double ODESolver::integrate(std::function<double(int32_t j, double z)> &integrand, int32_t j) const
{
    auto in = (int32_t) m_context.in;
    double res = 0;

    for (int32_t i = 0; i < in; ++i)
    {
        res += integrand(j, (double) (i + 1) / in) / in;
    }

    return res;
}

ODESolver::ODESolver(const Context &context)
{
    m_context = context;
    m_newton_solver = NewtonSolver(m_context);
    m_utils = Utils(m_context);

    dim = context.x_0.size();
}

void ODESolver::ComputeInitialAdamsX(double &t, std::vector<double> &x)
{
    for (uint32_t j = 0; j < m_context.adams_order; ++j)
    {

        // TODO optimize this
        std::vector<double> k_1 = k1(t, x);
        std::vector<double> k_2 = k2(t, x, k_1);
        std::vector<double> k_3 = k3(t, x, k_2);
        std::vector<double> k_4 = k4(t, x, k_3);

        for (uint32_t i = 0; i < x.size(); ++i)
        {
            x.at(i) += 1.0 / 6.0 * (k_1.at(i) + 2 * k_2.at(i) + 2 * k_3.at(i) + k_4.at(i));
        }

        Result.emplace_back(t, x);

        t += m_context.h;
    }
}

void ODESolver::rosenbrock()
{
    if (!m_context.f_autonomous)
    {
        std::cerr << "Autonomous ODE is not defined!" << std::endl;
        std::exit(1);
    }

    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    // Create Jacobi matrix
    std::vector<std::vector<double>> J(dim);
    for (auto &item: J)
    {
        item = std::vector<double>(dim);
    }

    // Magic numbers
    double alpha = 1.077;
    double betta = -0.372;
    double gamma = -0.577;

    while (t < m_context.t_end)
    {
        // Calculate Jacobi matrix
        m_utils.derive(J, m_context.f_autonomous, x);

        std::vector<std::vector<double>> JJ = algebra::multiply(J, J);

        // Left part
        algebra::multiply(JJ, -betta * m_context.h * m_context.h);
        algebra::multiply(J, -alpha * m_context.h);
        algebra::summ(JJ, J);
        algebra::summ(JJ, algebra::one(dim));
        algebra::multiply(JJ, 1.0 / m_context.h);

        //Right part
        std::vector<double> jj = m_context.f_autonomous(
                algebra::summ(x, algebra::multiply(m_context.f_autonomous(x), gamma * m_context.h))
                );

        // Construct system and solve

        for (algebra::size_type i = 0; i < dim; ++i)
        {
            JJ.at(i).push_back(jj.at(i));
        }

        std::vector<double> tmp(x);

        m_utils.solve(JJ, x);

        for (algebra::size_type i = 0; i < dim; ++i)
        {
            x.at(i) += tmp.at(i);
        }

        Result.emplace_back(t, x);

        t += m_context.h;
    }
}

void ODESolver::precor()
{
    if (!m_context.f)
    {
        std::cerr << "ODE is not defined!" << std::endl;
        std::exit(1);
    }

    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    computeA();
    computeB();

    // Compute n initial values of x
    ComputeInitialAdamsX(t, x);

    // Create n initial values of f
    std::vector<std::vector<double>> fs = ComputeInitialAdamsF();

    // Function to solve using Newton's method
    std::function<std::vector<double>(std::vector<double>)> f =
            [&fs, &t, this, &x](std::vector<double> y) -> std::vector<double>
            {
                std::vector<double> res = m_context.f(t, y);

                for (double &re: res)
                {
                    re *= B[0];
                }

                for (algebra::size_type j = 0; j < m_context.adams_order - 1; ++j)
                {
                    for (algebra::size_type i = 0; i < dim; ++i)
                    {
                        res[i] += B[m_context.adams_order - 1 - j] * fs[j + 1][i];
                    }
                }

                for (double &re: res)
                {
                    re *= m_context.h;
                }

                for (algebra::size_type i = 0; i < dim; ++i)
                {
                    res[i] += x[i] - y[i];
                }

                return res;
            };

    while (t < m_context.t_end)
    {
        // Create predictor using ae

        std::vector<double> guess = x;

        for (algebra::size_type j = 0; j < m_context.adams_order; ++j)
        {
            for (algebra::size_type i = 0; i < dim; ++i)
            {
                guess.at(i) += A.at(m_context.adams_order - 1 - j) * fs.at(j).at(i) * m_context.h;
            }
        }

        // Correct using ei
        x = m_newton_solver.solve_newton(f, guess);

        // Add iteration to the result and increment t
        Result.emplace_back(t, x);

        t += m_context.h;

        // Shift fs and compute next
        for (algebra::size_type i = 0; i < fs.size() - 1; ++i)
        {
            fs.at(i) = std::move(fs.at(i + 1));
        }
        fs.at(fs.size() - 1) = m_context.f(t, x);
    }
}

std::vector<std::vector<double>> ODESolver::ComputeInitialAdamsF()
{
    std::vector<std::vector<double>> fs(m_context.adams_order);

    for (uint32_t i = 0; i < Result.size(); ++i)
    {
        uint32_t index = Result.size() - 1 - i;
        fs.at(i) = m_context.f(
                std::get<0>(Result.at(index)),
                std::get<1>(Result.at(index))
        );
    }

    return fs;
}

