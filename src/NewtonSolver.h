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

namespace newton{
    typedef std::function<std::vector<double>(std::vector<double>)> function;

    /**
     * Returns
     */
    double residual(const vector &a, const vector &b);

    /**
     * Computes derivative (Jordan matrix)
     * @param m result Jordan matrix
     * @param f function to derive
     * @param x point in witch to compute derivative
     */
    void derive(matrix &m, const function& f, const vector &x);

    /**
     * Solves linear equation system
     * @param mx matrix with left and right parts of the system
     * @param res vector of results
     */
    void solve(matrix &mx, vector &res);

    /**
     * Solves system of non-linear equations f(x) = 0 using Newton's method
     * @return vector of results
     */
    vector solve_newton(const function& f,
                        const vector &initial_guess, uint32_t max_iterations);
}

#endif //NUMERICAL_TASK_9_NEWTONSOLVER_H
