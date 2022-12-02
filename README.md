# ODE Solver

# Features

Can solve systems of ordinary differential equation by using methods:

* Runge-Kutta fourth order method
* Adams extrapolation n-th order
* Adams interpolation n-th order
* Rosenbrock method for autonomous systems
* Predictor-Corrector method using Adams extrapolation and interpolation methods

# How to use

Create ODE using `std::function` and create context:

```c++
std::function<std::vector<double>(double t, std::vector<double> x)> f = 
        [](double t, std::vector<double> x) -> std::vector<double>{
    std::vector<double> res(x.size());

    /*
     * Your system. For example:
     */
    
    res[0] = 10.0 * (x[1] - x[0]);
    res[1] = x[0] * (28.0 - x[2]) - x[1];
    res[2] = x[0] * x[1] - 8.0 * x[2] / 3.0;

    return res;
};



Context context = Context(f)
```

Then set the initial values, start and begin of integration and other parameters.

```c++
context.x_0 = {10, 10, 10};
context.t_begin = 0.0;
context.t_end = 100.0;
context.h = 5e-5; // Integration step
context.adams_order = 10; // Order of Adams methods
```

Then you can create instance of `ODESolver`:

```c++
ODESolver odeSolver = ODESolver(context);
```

Run integration by one of the methods: `odeSolver.rk()`, `odeSolver.ae()`, `odeSolver.ai()`. Then retrieve results from `odeSolver.Results`:

```c++
odeSolver.ae();

std::ofstream ae("ae.dat");

for (auto &step: odeSolver.Result)
{
    auto [ t, x ] = step;

    ae << std::fixed << std::setprecision(10) << t << " ";

    for (auto &item: x)
    {
        ae << item << " ";
    }

    ae << std::endl;
}

ae.close();
```

`odeSolver.Result` has type of `std::vector<std::tuple<double, std::vector<double>>>`.
