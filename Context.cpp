//
// Created by Severin on 06.11.2022.
//

#include "Context.h"

/**
 * Vector of initial values
 */
std::vector<double> Context::x_0 = {0, 0};
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
uint32_t Context::n = 10;

std::vector<double> Context::f(double t, std::vector<double> x)
{
    for(auto & item : x){
        item *= item;
        item += t;
    }
    
    return x;
}
