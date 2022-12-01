#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <tuple>

#include "ODESolver.h"

void test_rk(const Context &context);

void test_ae(const Context &context);

void test_ai(const Context &context);

void test_rosen(const Context &context);

void test_precor(const Context &context);

void out(std::vector<std::tuple<double, std::vector<double>>> &res, const std::string& filename);

int main()
{
    // Creating function, Context with function and creating ODESolver with Context
    std::function<std::vector<double>(double t, std::vector<double> x)> f =
            [](double t, std::vector<double> x) -> std::vector<double>
            {
                std::vector<double> res(x.size());

                double kappa = 10;

                res[0] = x[1];
                res[1] = kappa * (1 - x[0] * x[0]) * x[1] - x[0];

                return res;
            };

    std::function<std::vector<double>(std::vector<double> x)> f_autonomous = [](std::vector<double> x) -> std::vector<double>
    {
        std::vector<double> res(x.size());

        double kappa = 10;

        res[0] = x[1];
        res[1] = kappa * (1 - x[0] * x[0]) * x[1] - x[0];

        return res;
    };

    Context context = Context(f);
    context.f_autonomous = f_autonomous;
    context.x_0 = {1, 1};
    context.adams_order = 5;
    context.h = 1e-4;
    context.t_begin = 0;
    context.t_end = 10;

    std::thread rk(test_rk, context);
    std::thread ae(test_ae, context);
    std::thread ai(test_ai, context);
    std::thread rosen(test_rosen, context);
    std::thread precor(test_precor, context);

    rk.join();
    ae.join();
    ai.join();
    rosen.join();
    precor.join();

    return 0;
}

void out(std::vector<std::tuple<double, std::vector<double>>> &res, const std::string& filename){
    std::ofstream file(filename);

    for (auto &step: res)
    {

        auto [t, x] = step;

        file << std::fixed << std::setprecision(20) << t << " ";

        for (auto &item: x)
        {
            file << item << " ";
        }

        file << std::endl;
    }

    file.close();
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

void test_rosen(const Context &context){
    ODESolver odeSolver = ODESolver(context);

    odeSolver.rosenbrock();

    out(odeSolver.Result, "rosen.dat");
}

void test_precor(const Context &context){
    ODESolver odeSolver = ODESolver(context);

    odeSolver.pc();

    out(odeSolver.Result, "precor.dat");
}
