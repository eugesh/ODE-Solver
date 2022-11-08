#include <iostream>
#include <fstream>
#include <iomanip>

#include "Integrator.h"

int main()
{
    Integrator integrator = Integrator();

    integrator.rk();

    std::ofstream rk("rk.dat");

    for (auto &step: integrator.Result)
    {
        auto [ t, x ] = step;

        rk << std::fixed << std::setprecision(10) << t << " ";

        for (auto &item: x)
        {
            rk << item << " ";
        }

        rk << std::endl;
    }

    rk.close();

    return 0;
}
