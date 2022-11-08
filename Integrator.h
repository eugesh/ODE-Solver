//
// Created by Severin on 08.11.2022.
//

#ifndef NUMERICAL_TASK_9_INTEGRATOR_H
#define NUMERICAL_TASK_9_INTEGRATOR_H

#include <vector>
#include <tuple>

#include "Context.h"

class Integrator
{
public:
    void rk();
    void ai();
    void ae();

    std::vector<std::tuple<double, std::vector<double>>> Result;

private:
    static std::vector<double> k1(double t, std::vector<double> x);
    static std::vector<double> k2(double t, std::vector<double> x, std::vector<double> k_1);
    static std::vector<double> k3(double t, std::vector<double> x, std::vector<double> k_2);
    static std::vector<double> k4(double t, std::vector<double> x, std::vector<double> k_3);
};


#endif //NUMERICAL_TASK_9_INTEGRATOR_H
