//
// Created by Severin on 06.11.2022.
//

#ifndef NUMERICAL_TASK_9_CONTEXT_H
#define NUMERICAL_TASK_9_CONTEXT_H

#include <vector>
#include <cstdint>
#include <functional>
#include <cmath>

class Context
{
public:
    explicit Context(std::function<std::vector<double>(double t, std::vector<double> x)> &f);

    explicit Context(std::function<std::vector<double>(std::vector<double> x)> &f_autonomous);

    Context() = default;

    /**
     * Right part of ODE x' = f(t, x)
     * @return x'
     */
    std::function<std::vector<double>(double t, std::vector<double> x)> f;
    /**
     * Right part of autonomous ODE x' = f(x). Used in Rosenbrock method for solving stiff equations.
     * @return x'
     */
    std::function<std::vector<double>(std::vector<double> x)> f_autonomous;
    /**
     * Vector of initial values. For system of n equations n initial values
     */
    std::vector<double> x_0 = {-10};
    /**
     * RK integration step
     */
    double h = 10e-5;
    /**
     * Begin of an integrating
     */
    double t_begin = 0.0;
    /**
     * End of a integrating
     */
    double t_end = 1.0;
    /**
     * Order of Adams interpolation / extrapolation method
     */
    uint32_t n = 4;
    /**
     * Number of integration steps in calculating Anj and Bnj for Adams methods
     */
    uint32_t in = 100000;

    double derive_step = 10e-8;
    double newton_precision = 10e-8;
    uint32_t newton_max_iterations = 1000;
};

#endif //NUMERICAL_TASK_9_CONTEXT_H
