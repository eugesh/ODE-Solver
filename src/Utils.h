//
// Created by Severin on 24.11.2022.
//

#ifndef NUMERICAL_TASK_9_MATHUTILS_H
#define NUMERICAL_TASK_9_MATHUTILS_H

#include "Context.h"
#include "Algebra.h"

class Utils
{
public:
    explicit Utils(const Context &context);
    Utils() = delete;

    /**
     * Computes derivative (Jordan matrix)
     * @param m result Jordan matrix
     * @param f function to derive
     * @param x point in witch to compute derivative
     */
    void derive(std::vector<std::vector<double>> &mx, std::function<std::vector<double>(std::vector<double>)> &f,
                const std::vector<double> &x) const;

    /**
     * Solves linear equation system
     * @param mx matrix with left and right parts of the system
     * @param res vector of results
     */
    void solve(std::vector<std::vector<double>> &mx, std::vector<double> &res);

    double residual(const std::vector<double> &a, const std::vector<double> &b);

private:
    Context m_context;
};

#endif //NUMERICAL_TASK_9_MATHUTILS_H
