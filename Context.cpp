//
// Created by Severin on 06.11.2022.
//

#include "Context.h"

#include <cmath>

/**
 * Vector of initial values
 */
std::vector<double> Context::x_0 = {-3, -1};
/**
 * RK integration step
 */
double Context::h = 10e-5;
/**
 * End of a integrating interval
 */
double Context::t_end = 3.0;
/**
 * Order of Adams interpolation / extrapolation method
 */
uint32_t Context::n = 10;

std::vector<double> Context::f(double t, std::vector<double> x)
{
    for(auto & item : x){
        item = pow(item, 2);
    }
    
    return x;
}
