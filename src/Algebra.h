//
// Created by Severin on 18.11.2022.
//

#ifndef NUMERICAL_TASK_9_ALGEBRA_H
#define NUMERICAL_TASK_9_ALGEBRA_H

#include <vector>
#include <functional>

namespace algebra
{
    typedef std::vector<double>::size_type size_type;

    /**
     * Vector difference
     */
    std::vector<double> subtract(const std::vector<double> &a, const std::vector<double> &b);

    /**
     * Matrix by vector multiplication
     */
    std::vector<double> multiply(const std::vector<std::vector<double>> &m, const std::vector<double> &x);

    /**
     * Vector by number division
     */
    std::vector<double> divide(const std::vector<double> &a, double value);

    /**
     * Creates n x n + 1 matrix from n x n matrix and n x 1 vector
     * @param m initial and resulted matrix
     * @param x vector
     */
    void plug_vector(std::vector<std::vector<double>> &m, const std::vector<double> &x);
}

#endif //NUMERICAL_TASK_9_ALGEBRA_H
