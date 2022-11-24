//
// Created by Severin on 24.11.2022.
//

#include "Utils.h"

Utils::Utils(Context &context)
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
        l[column] += m_context.newton_derive_step;
        r[column] -= m_context.newton_derive_step;

        std::vector<double> d = algebra::divide(algebra::subtract(f(l), f(r)), 2 * m_context.newton_derive_step);

        for (algebra::size_type row = 0; row < n; ++row)
        {
            mx[row][column] = d[row];
        }

        l[column] -= m_context.newton_derive_step;
        r[column] += m_context.newton_derive_step;
    }
}

Utils::Utils()
{
    m_context = Context();
}

