//
// Created by Severin on 08.11.2022.
//

#ifndef NUMERICAL_TASK_9_ODESOLVER_H
#define NUMERICAL_TASK_9_ODESOLVER_H

#include <vector>
#include <tuple>

#include "Context.h"

class ODESolver
{
public:
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

    std::vector<std::tuple<double, std::vector<double>>> Result;

private:
    /**
     * Methods for calculating KR4 coefficients
     */
    static std::vector<double> k1(double t, std::vector<double> x);
    static std::vector<double> k2(double t, std::vector<double> x, std::vector<double> k_1);
    static std::vector<double> k3(double t, std::vector<double> x, std::vector<double> k_2);
    static std::vector<double> k4(double t, std::vector<double> x, std::vector<double> k_3);

    /**
     * Methods for computing A and B
     */

    void computeA();
    void computeB();

    /**
     * Coefficients Anj for Adams extrapolation method
     */
    std::vector<double> A = std::vector<double>(Context::n);

    /**
     * Coefficients Bnj for Adams interpolation method
     */
    std::vector<double> B = std::vector<double>(Context::n);

    /**
     * Factorial
     */
    static int32_t factorial(int32_t x);

    /**
     * Integration method from 0 to 1 with 1000 intervals
     */

    static double integrate(double (*integrand)(int32_t j, double z), int32_t j);

    /**
     * Integrand for Adams method (ea)
     */
    static double integrandForA(int32_t j, double z);

    /**
     * Integrand for Adams method (ea)
     */
    static double integrandForB(int32_t j, double z);
};


#endif //NUMERICAL_TASK_9_ODESOLVER_H
