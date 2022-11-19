//
// Created by Severin on 08.11.2022.
//

#include "ODESolver.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <functional>

#include "NewtonSolver.h"

void ODESolver::rk()
{
    Result.clear();
    double t = Context::x_0.at(0); // Get t_0 from x_0
    std::vector<double> x = Context::x_0;

    x.erase(x.begin()); // Remove t_0 from x_0

    while (t < Context::t_end)
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

        t += Context::h;
    }
}

std::vector<double> ODESolver::k1(double t, std::vector<double> x)
{
    std::vector<double> res = Context::f(t, std::move(x));

    for (auto &item: res)
    {
        item *= Context::h;
    }

    return res;
}

std::vector<double> ODESolver::k2(double t, std::vector<double> x, std::vector<double> k_1)
{
    for (uint32_t i = 0; i < x.size(); ++i)
    {
        x.at(i) += k_1.at(i) / 2;
    }

    std::vector<double> res = Context::f(t + Context::h / 2, std::move(x));

    for (auto &item: res)
    {
        item *= Context::h;
    }

    return res;

}

std::vector<double> ODESolver::k3(double t, std::vector<double> x, std::vector<double> k_2)
{
    for (uint32_t i = 0; i < x.size(); ++i)
    {
        x.at(i) += k_2.at(i) / 2;
    }

    std::vector<double> res = Context::f(t + Context::h / 2, std::move(x));

    for (auto &item: res)
    {
        item *= Context::h;
    }

    return res;

}

std::vector<double> ODESolver::k4(double t, std::vector<double> x, std::vector<double> k_3)
{
    for (uint32_t i = 0; i < x.size(); ++i)
    {
        x.at(i) += k_3.at(i);
    }

    std::vector<double> res = Context::f(t + Context::h, std::move(x));

    for (auto &item: res)
    {
        item *= Context::h;
    }

    return res;

}

void ODESolver::ae()
{
    Result.clear();
    double t = Context::x_0.at(0); // Get t_0 from x_0
    std::vector<double> x = Context::x_0;

    x.erase(x.begin()); // Remove t_0 from x_0

    computeA();

    // Compute n initial values of x

    for (uint32_t j = 0; j < Context::n; ++j)
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

        t += Context::h;
    }

    // Run Adams method

    // Create n initial values of f
    std::vector<std::vector<double>> fs(Context::n);
    for (uint32_t i = 0; i < Context::n; ++i)
    {
        uint32_t index = Result.size() - 1 - i;
        fs.at(i) = Context::f(
                std::get<0>(Result.at(index)),
                std::get<1>(Result.at(index))
        );
    }

    while (t < Context::t_end)
    {

        for (uint32_t j = 0; j < Context::n; ++j)
        {
            for (uint32_t i = 0; i < x.size(); ++i)
            {
                x.at(i) += A.at(j) * fs.at(j).at(i) * Context::h;
            }
        }

        Result.emplace_back(t, x);

        t += Context::h;

        for (uint32_t i = 0; i < fs.size() - 1; ++i)
        {
            fs.at(i) = std::move(fs.at(i + 1));
        }

        fs.at(fs.size() - 1) = Context::f(t, x);
    }
}

void ODESolver::ai()
{
    Result.clear();
    double t = Context::x_0.at(0); // Get t_0 from x_0
    std::vector<double> x = Context::x_0;

    x.erase(x.begin()); // Remove t_0 from x_0

    computeB();

    // Compute n initial values of x

    for (uint32_t j = 0; j < Context::n - 1; ++j)
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

        t += Context::h;
    }

    // Run Adams method

    // Create n initial values of f
    std::vector<std::vector<double>> fs(Context::n - 1);
    for (uint32_t i = 0; i < Context::n - 1; ++i)
    {
        uint32_t index = Result.size() - 1 - i;
        fs.at(i) = Context::f(
                std::get<0>(Result.at(index)),
                std::get<1>(Result.at(index))
        );
    }

    while (t < Context::t_end)
    {

        std::function<std::vector<double>(std::vector<double>)> f =
                [&fs, t, this](std::vector<double> y) -> std::vector<double>
                {
                    std::vector<double> res = Context::f(t, y);

                    for (double &re: res)
                    {
                        re *= B.at(Context::n - 1);
                    }

                    for (uint32_t j = 0; j < Context::n - 1; ++j)
                    {
                        for (uint32_t i = 0; i < res.size(); ++i)
                        {
                            res.at(i) += B.at(j) * fs.at(j).at(i);
                        }
                    }

                    for (double &re: res)
                    {
                        re *= Context::h;
                    }

                    for (uint32_t i = 0; i < res.size(); ++i)
                    {
                        res.at(i) += fs.at(fs.size() - 1).at(i) - y.at(i);
                    }

                    return res;
                };

        x = newton::solve_newton(f, x, 1000);
        Result.emplace_back(t, x);

        t += Context::h;

        for (uint32_t i = 0; i < fs.size() - 1; ++i)
        {
            fs.at(i) = std::move(fs.at(i + 1));
        }

        fs.at(fs.size() - 1) = Context::f(t, x);
    }
}

void ODESolver::computeA()
{
    //TODO optimize this
    for (int32_t j = 0; j < (int32_t) Context::n; ++j)
    {
        A.at(j) = pow(-1.0, j) / (factorial(j) * factorial((int32_t) Context::n - 1 - j));

        A.at(j) *= integrate(integrandForA, j);
    }
}

void ODESolver::computeB()
{
    //TODO optimize this
    for (int32_t j = -1; j < (int32_t) Context::n - 1; ++j)
    {
        A.at(j + 1) = pow(-1.0, j + 1) / (factorial(j + 1) * factorial((int32_t) Context::n - 2 - j));

        A.at(j + 1) *= integrate(integrandForB, j);
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

double ODESolver::integrate(double (*integrand)(int32_t, double), int32_t j)
{
    auto n = (int32_t) Context::in;
    double res = 0;

    for (int32_t i = 0; i < n; ++i)
    {
        res += integrand(j, (double) (i + 1) / n) / n;
    }

    return res;
}

double ODESolver::integrandForA(int32_t j, double z)
{
    double res = 1;

    for (int32_t i = 0; i < (int32_t) Context::n; ++i)
    {
        if (i == j)
        {
            continue;
        }
        res *= z + i;
    }

    return res;
}

double ODESolver::integrandForB(int32_t j, double z)
{
    double res = 1;

    for (int32_t i = -1; i < (int32_t) Context::n - 1; ++i)
    {
        if (i == j)
        {
            continue;
        }
        res *= z + i;
    }

    return res;
}
