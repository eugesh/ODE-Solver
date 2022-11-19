//
// Created by Severin on 06.11.2022.
//

#include "Context.h"

#include <cmath>

/**
 * Vector of initial values
 */
std::vector<double> Context::x_0 = {-1, -10};
/**
 * RK integration step
 */
double Context::h = 10e-5;
/**
 * End of a integrating interval
 */
double Context::t_end = 1.0;
/**
 * Order of Adams interpolation / extrapolation method
 */
uint32_t Context::n = 4;
/**
 * Number of integration steps in calculating Anj and Bnj for Adams methods
 */
uint32_t Context::in = 1000;

double Context::newton_derive_step = 10e-5;
double Context::newton_precision = 10e-5;

std::vector<double> Context::f(double t, std::vector<double> x)
{
    for(auto & item : x){
        item = 3.0 * pow(item * item, 0.333);
    }
    
    return x;
}
