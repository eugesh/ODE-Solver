#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>

#include "ODESolver.h"

void test_rk(const Context &context);

void test_ae(const Context &context);

void test_ai(const Context &context);

void out(std::vector<std::tuple<double, std::vector<double>>> &res, const std::string& filename);

int main()
{
    // Creating function, Context with function and creating ODESolver with Context
    std::function<std::vector<double>(double t, std::vector<double> x)> f =
            [](double t, std::vector<double> x) -> std::vector<double>
            {
                std::vector<double> res(x.size());

                res[0] = 10.0 * (x[1] - x[0]);
                res[1] = x[0] * (28.0 - x[2]) - x[1];
                res[2] = x[0] * x[1] - 8.0 * x[2] / 3.0;

                return res;
            };

    std::function<std::vector<double>(std::vector<double> x)> f_autonomous = [](std::vector<double> x) -> std::vector<double>
    {
        std::vector<double> res(x.size());

        res[0] = -0.05 * x[0] + 1e4 * x[1] * x[2];
        res[1] = 0.05 * x[0] - 1e4 * x[1] * x[2] - 1e7 * x[1];
        res[2] = 1e7 * x[1];

        return res;
    };

    Context context = Context(f);
    context.f_autonomous = f_autonomous;
    context.x_0 = {10, 10, 10};
    context.n = 10;
    context.h = 5e-5;
    context.t_end = 25;

    // Run calculation
    std::thread rk(test_rk, context);
    std::thread ae(test_ae, context);
    std::thread ai(test_ai, context);

    rk.join();
    ae.join();
    ai.join();

    return 0;
}

void out(std::vector<std::tuple<double, std::vector<double>>> &res, const std::string& filename){
    std::ofstream rk(filename);

    for (auto &step: res)
    {
        auto [t, x] = step;

        rk << std::fixed << std::setprecision(10) << t << " ";

        for (auto &item: x)
        {
            rk << item << " ";
        }

        rk << std::endl;
    }

    rk.close();
}

void test_rk(const Context &context)
{
    ODESolver odeSolver = ODESolver(context);

    odeSolver.rk();

    out(odeSolver.Result, "rk.dat");
}

void test_ae(const Context &context)
{
    ODESolver odeSolver = ODESolver(context);

    odeSolver.ae();

    out(odeSolver.Result, "ae.dat");
}

void test_ai(const Context &context)
{
    ODESolver odeSolver = ODESolver(context);

    odeSolver.ai();

    out(odeSolver.Result, "ai.dat");
}
