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
    std::vector<double> difference(const std::vector<double> &a, const std::vector<double> &b);

    /**
     * Vector by number multiplication
     */
    std::vector<double> multiply(const std::vector<double> &a, double value);


    /**
     * Vector summ
     */
    std::vector<double> summ(const std::vector<double> &a, const std::vector<double> &b);


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

    /**
     * Matrix summ
     */
    void summ(std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b);

    /**
     * Matrix difference
     */
    void difference(std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b);


    /**
     * Matrix by number multiplication
     */
    void multiply(std::vector<std::vector<double>> &a, double value);

    /**
     * Matrix by matrix multiplication
     */
    std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b);

    /**
     * Identity matrix
     */
    std::vector<std::vector<double>> one(size_type n);

    /**
     * Empty matrix
     */
    std::vector<std::vector<double>> zero(size_type n);

}

#endif //NUMERICAL_TASK_9_ALGEBRA_H
