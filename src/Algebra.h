//
// Created by Severin on 18.11.2022.
//

#ifndef NUMERICAL_TASK_9_ALGEBRA_H
#define NUMERICAL_TASK_9_ALGEBRA_H

#include <vector>
#include <functional>

namespace algebra{
    typedef std::vector<double> vector;
    typedef std::vector<std::vector<double>> matrix;
    typedef vector::size_type size_type;
    typedef std::vector<size_type> vector_int;

    /**
     * Vector difference
     */
    vector subtract(const vector &a, const vector &b);

    /**
     * Matrix by vector multiplication
     */
    vector multiply(const matrix &m, const vector &x);

    /**
     * Vector by number division
     */
    vector divide(const vector &a, double value);

    /**
     * Creates matrix
     * @param m result matrix
     * @param height matrix height
     * @param with matrix with
     */
    void create(matrix &m, size_type height, size_type with);

    /**
     * Creates n x n + 1 matrix from n x n matrix and n x 1 vector
     * @param m initial and resulted matrix
     * @param x vector
     */
    void plug_vector(matrix &m, const vector &x);
}

#endif //NUMERICAL_TASK_9_ALGEBRA_H
