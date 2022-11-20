#include <iostream>
#include <fstream>
#include <iomanip>

#include "ODESolver.h"

int main()
{
    // Creating function, Context with function and creating ODESolver with Context
    std::function<std::vector<double>(double t, std::vector<double> x)> f =
            [](double t, std::vector<double> x) -> std::vector<double>{
        std::vector<double> res(x.size());

        res.at(0) = sin(x.at(0) / x.at(1));
        res.at(1) = x.at(0) / x.at(2) + t;
        res.at(2) = exp(x.at(1) - x.at(2));

        return res;
    };

    Context context = Context(f);
    context.x_0 = {-6, 1, 5, 5};
    context.n = 5;
    context.h = 10e-5;
    context.t_end = 6;
    ODESolver odeSolver = ODESolver(context);

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
