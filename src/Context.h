//
// Created by Severin on 06.11.2022.
//

#ifndef NUMERICAL_TASK_9_CONTEXT_H
#define NUMERICAL_TASK_9_CONTEXT_H

#include <vector>
#include <cstdint>
#include <functional>
#include <cmath>

class Context{
public:
    explicit Context(std::function<std::vector<double>(double t, std::vector<double> x)> &f);
    Context();
    /**
     * Right part of ODE x' = f(t, x)
     * @return x'
     */
    std::function<std::vector<double>(double t, std::vector<double> x)> f;
    /**
     * Vector of initial values
     */
    std::vector<double> x_0 = {-1, -10};
    /**
     * RK integration step
     */
    double h = 10e-5;
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
    uint32_t in = 1000;

    double newton_derive_step = 10e-5;
    double newton_precision = 10e-5;
};

#endif //NUMERICAL_TASK_9_CONTEXT_H
