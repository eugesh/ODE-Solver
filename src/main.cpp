#include <iostream>
#include <fstream>
#include <iomanip>

#include "ODESolver.h"

int main()
{
    ODESolver odeSolver = ODESolver();

    // rk

    odeSolver.rk();

    std::ofstream rk("rk.dat");

    for (auto &step: odeSolver.Result)
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

    // ae

    odeSolver.ae();

    std::ofstream ae("ae.dat");

    for (auto &step: odeSolver.Result)
    {
        auto [ t, x ] = step;

        ae << std::fixed << std::setprecision(10) << t << " ";

        for (auto &item: x)
        {
            ae << item << " ";
        }

        ae << std::endl;
    }

    ae.close();

    // ai

    odeSolver.ai();

    std::ofstream ai("ai.dat");

    for (auto &step: odeSolver.Result)
    {
        auto [ t, x ] = step;

        ai << std::fixed << std::setprecision(10) << t << " ";

        for (auto &item: x)
        {
            ai << item << " ";
        }

        ai << std::endl;
    }

    ai.close();

    return 0;
}
