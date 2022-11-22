//
// Created by Severin on 18.11.2022.
//

#ifndef NUMERICAL_TASK_9_NEWTONSOLVER_H
#define NUMERICAL_TASK_9_NEWTONSOLVER_H

#include <vector>
#include <cstdint>

#include "Algebra.h"
#include "Context.h"

using namespace algebra;

class NewtonSolver{
public:
    explicit NewtonSolver(const Context& context);
    NewtonSolver();

    /**
     * Solves system of non-linear equations f(x) = 0 using Newton's method
     * @return vector of results
     */
    std::vector<double> solve_newton(std::function<std::vector<double>(std::vector<double>)>& f,
                        const std::vector<double> &initial_guess, uint32_t max_iterations);
private:
    Context m_context;

    static double residual(const std::vector<double> &a, const std::vector<double> &b);

    /**
     * Computes derivative (Jordan matrix)
     * @param m result Jordan matrix
     * @param f function to derive
     * @param x point in witch to compute derivative
     */
    void derive(std::function<std::vector<double>(std::vector<double>)>& f, const std::vector<double> &x);

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
