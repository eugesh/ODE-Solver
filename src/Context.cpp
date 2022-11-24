//
// Created by Severin on 06.11.2022.
//

#include "Context.h"

Context::Context()
{
    f = [](double t, std::vector<double> x) -> std::vector<double>
    {
        for (auto &item: x)
        {
            item = 3.0 * pow(item * item, 0.333);
        }

        return x;
    };

    f_autonomous = [](std::vector<double> x) -> std::vector<double>
    {
        std::vector<double> res(x.size());

        res[0] = -0.05 * x[0] + 1e4 * x[1] * x[2];
        res[1] = 0.05 * x[0] - 1e4 * x[1] * x[2] - 1e7 * x[1];
        res[2] = 1e7 * x[1];

        return res;
    };
}

Context::Context(std::function<std::vector<double>(double, std::vector<double>)> &t_f)
{
    f = t_f;
}

Context::Context(std::function<std::vector<double>(std::vector<double>)> &t_f_autonomous)
{
    f_autonomous = t_f_autonomous;
}
