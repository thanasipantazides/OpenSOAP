# Numerical Methods
Some overview of numerical methods used in this package.

## Integration
For orbital and attitude dynamics propagation, I rely on the [RK4 numerical integration scheme](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods). This is not the "best" way to integrate ordinary differential equations ([this package](https://github.com/SciML/DifferentialEquations.jl) will tell you the best way), but it is really simple to implement, and allows me to get into the integrator directly and mess with states. 

The ability to directly mess with the state vector after each integration step is handy for propagating attitude dynamics.