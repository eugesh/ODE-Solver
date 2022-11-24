//
// Created by Severin on 08.11.2022.
//

#ifndef NUMERICAL_TASK_9_ODESOLVER_H
#define NUMERICAL_TASK_9_ODESOLVER_H

#include <vector>
#include <tuple>

#include "Context.h"
#include "NewtonSolver.h"

class ODESolver
{
public:
    explicit ODESolver(const Context &context);

    /**
     * Runge–Kutta 4 integrator
     */
    void rk();

    /**
     * Adams extrapolation integrator
     */
    void ae();

    /**
     * Adams interpolation integrator
     */
    void ai();

    /**
     * Rosenbrock interrogator
     */
    void rosenbrock();

    /**
     * Predictor–corrector interrogator
     */
    void pc();

    std::vector<std::tuple<double, std::vector<double>>> Result;

private:
    Context m_context;
    NewtonSolver m_newton_solver;

    /**
     * Methods for calculating KR4 coefficients
     */
    [[nodiscard]] std::vector<double> k1(double t, std::vector<double> x) const;

    [[nodiscard]] std::vector<double> k2(double t, std::vector<double> x, std::vector<double> k_1) const;

    [[nodiscard]] std::vector<double> k3(double t, std::vector<double> x, std::vector<double> k_2) const;

    [[nodiscard]] std::vector<double> k4(double t, std::vector<double> x, std::vector<double> k_3) const;

    /**
     * Methods for computing A and B
     */

    void computeA();

    void computeB();

    /**
     * Compute initial values for adams methods extrapolating methods
     */
    void ComputeInitialAE(double &t, std::vector<double> &x);

    /**
     * Compute initial values for adams methods interpolating methods
     */
    void ComputeInitialAI(double &t, std::vector<double> &x);

    /**
     * Coefficients Anj for Adams extrapolation method
     */
    std::vector<double> A;

    /**
     * Coefficients Bnj for Adams interpolation method
     */
    std::vector<double> B;

    /**
     * Factorial
     */
    static int32_t factorial(int32_t x);

    /**
     * Integration method from 0 to 1 with 1000 intervals
     */

    double integrate(std::function<double(int32_t j, double z)> &integrand, int32_t j) const;

    /**
     * Integrand for Adams method (ea)
     */
    std::function<double(int32_t j, double z)> integrandForA = [this](int32_t j, double z) -> double
    {
        double res = 1;

        for (int32_t i = 0; i < (int32_t) m_context.n; ++i)
        {
            if (i == j)
            {
                continue;
            }
            res *= z + i;
        }

        return res;
    };

    /**
     * Integrand for Adams method (ia)
     */
    std::function<double(int32_t j, double z)> integrandForB = [this](int32_t j, double z) -> double
    {
        double res = 1;

        for (int32_t i = -1; i < (int32_t) m_context.n - 1; ++i)
        {
            if (i == j)
            {
                continue;
            }
            res *= z + i;
        }

        return res;
    };
};


#endif //NUMERICAL_TASK_9_ODESOLVER_H
