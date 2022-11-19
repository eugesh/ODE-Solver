//
// Created by Severin on 18.11.2022.
//

#include "Algebra.h"

namespace algebra{
    void create(matrix &m, size_type height, size_type with)
    {
        m = matrix(height);

        for (auto &item: m)
        {

            item = vector(with);
        }
    }

    vector multiply(const matrix &m, const vector &x)
    {
        size_type n = m.size();
        vector result = vector(n);

        for (size_type i = 0; i < n; ++i)
        {
            for (size_type j = 0; j < n; ++j)
            {
                result[i] += m[i][j] * x[j];
            }
        }

        return result;
    }

    vector subtract(const vector &a, const vector &b)
    {
        size_type n = a.size();
        vector result = vector(n);

        for (size_type i = 0; i < n; ++i)
        {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    vector divide(const vector &a, double value)
    {
        size_type n = a.size();
        vector result = vector(n);

        for (size_type i = 0; i < n; ++i)
        {
            result[i] = a[i] / value;
        }

        return result;
    }

    void plug_vector(matrix &m, const vector &x)
    {
        size_type n = m.size();

        for (size_type i = 0; i < n; ++i)
        {
            m[i][n] = x[i];
        }
    }
}
