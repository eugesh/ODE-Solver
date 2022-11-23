//
// Created by Severin on 08.11.2022.
//

#include "ODESolver.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <utility>
#include <iostream>

void ODESolver::rk()
{
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

        for (uint32_t i = 0; i < x.size(); ++i)
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
    for (uint32_t i = 0; i < x.size(); ++i)
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
    for (uint32_t i = 0; i < x.size(); ++i)
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
    for (uint32_t i = 0; i < x.size(); ++i)
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
    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    computeA();

    // Compute n initial values of x

    for (uint32_t j = 0; j < m_context.n; ++j)
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

    // Run Adams method

    // Create n initial values of f
    std::vector<std::vector<double>> fs(m_context.n);
    for (uint32_t i = 0; i < m_context.n; ++i)
    {
        uint32_t index = Result.size() - 1 - i;
        fs.at(i) = m_context.f(
                std::get<0>(Result.at(index)),
                std::get<1>(Result.at(index))
        );
    }

    while (t < m_context.t_end)
    {

        for (uint32_t j = 0; j < m_context.n; ++j)
        {
            for (uint32_t i = 0; i < x.size(); ++i)
            {
                x.at(i) += A.at(m_context.n - 1 -j) * fs.at(j).at(i) * m_context.h;
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
    Result.clear();

    double t = m_context.t_begin;
    std::vector<double> x = m_context.x_0;

    computeB();

    // Compute n initial values of x

    for (uint32_t j = 0; j < m_context.n - 1; ++j)
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

    // Run Adams method

    // Create n initial values of f
    std::vector<std::vector<double>> fs(m_context.n - 1);
    for (uint32_t i = 0; i < Result.size(); ++i)
    {
        uint32_t index = Result.size() - 1 - i;
        fs.at(i) = m_context.f(
                std::get<0>(Result.at(index)),
                std::get<1>(Result.at(index))
        );
    }

    std::function<std::vector<double>(std::vector<double>)> f =
            [&fs, &t, this, &x](std::vector<double> y) -> std::vector<double>
            {
                std::vector<double> res = m_context.f(t, y);

                for (double &re: res)
                {
                    re *= B[0];
                }

                for (uint32_t j = 0; j < fs.size(); ++j)
                {
                    for (uint32_t i = 0; i < res.size(); ++i)
                    {
                        res[i] += B[m_context.n - 1 - j] * fs[j][i];
                    }
                }

                for (double &re: res)
                {
                    re *= m_context.h;
                }

                for (uint32_t i = 0; i < res.size(); ++i)
                {
                    res[i] += x[i] - y[i];
                }

                return res;
            };

    while (t < m_context.t_end)
    {
        x = m_newton_solver.solve_newton(f, x, m_context.newton_max_iterations);
        Result.emplace_back(t, x);

        t += m_context.h;

        for (uint32_t i = 0; i < fs.size() - 1; ++i)
        {
            fs.at(i) = std::move(fs.at(i + 1));
        }

        fs.at(fs.size() - 1) = m_context.f(t, x);
    }
}

void ODESolver::computeA()
{
    A = std::vector<double>(m_context.n);

    //TODO optimize this
    for (int32_t j = 0; j < (int32_t) m_context.n; ++j)
    {
        A.at(j) = pow(-1.0, j) / (factorial(j) * factorial((int32_t) m_context.n - 1 - j));

        A.at(j) *= integrate(integrandForA, j);
    }
}

void ODESolver::computeB()
{
    B = std::vector<double>(m_context.n);

    //TODO optimize this
    for (int32_t j = -1; j < (int32_t) m_context.n - 1; ++j)
    {
        B.at(j + 1) = pow(-1.0, j + 1) / (factorial(j + 1) * factorial((int32_t) m_context.n - 2 - j));

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

double ODESolver::integrate(std::function<double(int32_t j, double z)>& integrand, int32_t j) const
{
    auto n = (int32_t) m_context.in;
    double res = 0;

    for (int32_t i = 0; i < n; ++i)
    {
        res += integrand(j, (double) (i + 1) / n) / n;
    }

    return res;
}

ODESolver::ODESolver(const Context& context)
{
    m_context = context;
    m_newton_solver = NewtonSolver(context);
}
