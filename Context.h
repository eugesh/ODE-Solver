//
// Created by Severin on 06.11.2022.
//

#ifndef NUMERICAL_TASK_9_CONTEXT_H
#define NUMERICAL_TASK_9_CONTEXT_H

#include <vector>
#include <cstdint>

class Context{
public:
    Context() = delete;
    /**
     * Right part of ODE x' = f(t, x)
     * @return x'
     */
    static std::vector<double> f(double t, std::vector<double> x);
    /**
     * Vector of initial values
     */
    static std::vector<double> x_0;
    /**
     * RK integration step
     */
    static double h;
    /**
     * End of a integrating
     */
    static double t_end;
    /**
     * Order of Adams interpolation / extrapolation method
     */
    static uint32_t n;
};

#endif //NUMERICAL_TASK_9_CONTEXT_H
