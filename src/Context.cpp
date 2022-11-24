//
// Created by Severin on 06.11.2022.
//

#include "Context.h"

Context::Context(std::function<std::vector<double>(double, std::vector<double>)> &t_f)
{
    f = t_f;
}

Context::Context(std::function<std::vector<double>(std::vector<double>)> &t_f_autonomous)
{
    f_autonomous = t_f_autonomous;
}
