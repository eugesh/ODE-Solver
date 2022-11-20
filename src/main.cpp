#include <iostream>
#include <fstream>
#include <iomanip>

#include "ODESolver.h"

int main()
{
    // Creating function, Context with function and creating ODESolver with Context
    std::function<std::vector<double>(double t, std::vector<double> x)> f =
            [](double t, std::vector<double> x) -> std::vector<double>{
        for(auto & item : x){
            item = 3.0 * pow(item * item, 0.333);
        }

        return x;
    };

    Context context = Context(f);
    context.n = 5;
    ODESolver odeSolver = ODESolver(Context());

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
