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
}

Context::Context(std::function<std::vector<double>(double, std::vector<double>)> &t_f)
{
    f = t_f;
}
