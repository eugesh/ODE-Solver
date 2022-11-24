//
// Created by Severin on 18.11.2022.
//

#ifndef NUMERICAL_TASK_9_NEWTONSOLVER_H
#define NUMERICAL_TASK_9_NEWTONSOLVER_H

#include <vector>
#include <cstdint>

#include "Algebra.h"
#include "Context.h"
#include "Utils.h"

using namespace algebra;

class NewtonSolver
{
public:
    explicit NewtonSolver(const Context &context);

    NewtonSolver() = delete;

    /**
     * Solves system of non-linear equations f(x) = 0 using Newton's method
     * @return vector of results
     */
    std::vector<double> solve_newton(std::function<std::vector<double>(std::vector<double>)> &f,
                                     const std::vector<double> &initial_guess, uint32_t max_iterations);

private:
    Context m_context;
    Utils m_utils = Utils(Context());

    static double residual(const std::vector<double> &a, const std::vector<double> &b);

    /**
     * Solves linear equation system
     * @param mx matrix with left and right parts of the system
     * @param res vector of results
     */
    static void solve(std::vector<std::vector<double>> &mx, std::vector<double> &res);

    /**
     * Matrix of derivatives n by n + 1: (derivatives | derivatives * previous - f(previous))
     */
    std::vector<std::vector<double>> m;
};

#endif //NUMERICAL_TASK_9_NEWTONSOLVER_H
