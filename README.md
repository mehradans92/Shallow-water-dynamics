# Shallow Water Dynamics
The [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations) (SWE), have been solved numerically using 4th order spatial and 3rd order temporal discretization in 1D and 2D. The adaptive timestep is chosen based on the [CFL criteria](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) throught the simulation. The model initiates with a gaussian bump in the middle which propegates with time into periodic boundaries under flow and no-flow conditions. The effect of topography (η) and earth's rotation (𝑓) have been included in the model. Higher 𝑓 results a stronger rotation force which can be observed in the velocity quiver plots as time increases (see animation). The key point here is that for strong Coriolis effects, variations in height and velocity are not radiated away to infinity, as in the non-rotating problem. Note that this rotation is anticlock-wise round the basin for positive values of 𝑓 in northern hemisphere and conversely in southern hemisphere.  
  
![](/Corriolis_effect.gif)
![](/topography_effect_2.gif)
![](/topography_effect.gif)

