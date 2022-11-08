//
// Created by Severin on 08.11.2022.
//

#include "Integrator.h"

void Integrator::rk()
{
    Result.clear();
    double t = Context::x_0.at(0); // Get t_0 from x_0
    std::vector<double> x = Context::x_0;

    x.erase(x.begin()); // Remove t_0 from x_0

    while (t < Context::t_end){

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

std::vector<double> Integrator::k1(double t, std::vector<double> x)
{
    std::vector<double> res = Context::f(t, std::move(x));

    for (auto &item: res)
    {
        item *= Context::h;
    }

    return res;
}

std::vector<double> Integrator::k2(double t, std::vector<double> x, std::vector<double> k_1)
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

std::vector<double> Integrator::k3(double t, std::vector<double> x, std::vector<double> k_2)
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

std::vector<double> Integrator::k4(double t, std::vector<double> x, std::vector<double> k_3)
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
